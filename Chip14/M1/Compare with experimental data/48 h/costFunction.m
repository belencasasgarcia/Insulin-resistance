function [error] = costFunction(parameters)


global modelName
global icOrig0
global OPTIONS
global COSTOPTIONS
global EXPDATA
global SIMTIME
global OPTVAR_INDEX
global param
global CUTOFF
global FID
global Optparam
global pNames
global parIndex

Optpar=Optparam;    % Create a copy of the Optparam variable since it will be modified                         
Optpar(OPTIONS.index_Optpar)=parameters;    % Update with current values of the parameters to be optimized   

% Calculates the model agreement with the data based on the size of the
% residuals.

%% Simulate

%% Simulation for the co-culture
try
    simData_islets_liver = feval(char(modelName), SIMTIME, icOrig0, Optpar, COSTOPTIONS);
    
catch simFail
    
    disp('Simulation crashed');
    error = inf;
    disp(param);
    disp(simFail)
    return;
end

%% Simulation for only liver
Optpar_liver(parIndex.i_V_islets)=0;

try
    simData_liver = feval(char(modelName), SIMTIME, icOrig0, Optpar_liver, COSTOPTIONS);
    
catch simFail
    
    disp('Simulation crashed');
    error = inf;
    disp(param);
    disp(simFail)
    return;
end

%% Simulation for only liver when insulin is added

% Simulate model with glucose 11 mM, added insulin 488.4311 mU/L
    
icOrig_hep_I=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets) ...
    488.4311*param(parIndex.i_V_m_hep) 488.4311*param(parIndex.i_V_m_islets)];
    
try
    simData_liver_insulin = feval(char(modelName), SIMTIME, icOrig_hep_I, Optpar_liver, COSTOPTIONS);
    
catch simFail
    
disp('Simulation crashed');
error = inf;
disp(param);
disp(simFail)

return;
end

tmpError = 0;

for n=OPTVAR_INDEX
    
% 1: Measured glucose with liver+islets (variable 1)
% 2: Measured insulin with liver+islets (variable 2)
% 3: Measured glucose with only liver   (variable 1, V_islets=0)
% 4: Measured insulin with only liver   (variable 2, V_islets=0)

if n==1 || n==2
    
    % Values from the co-culture simulation
    simData_values=simData_islets_liver.variablevalues(:,n);
    index_int=ismembertol(SIMTIME,EXPDATA.time{n},1e-4);

    tmpError = tmpError + sum(((simData_values(index_int)- EXPDATA.mean{n}).^2)./(EXPDATA.SD{n}).^2);

end

%% Glucose with only hepatocytes

if n==3% 3: Measured glucose (variable 1)

    simData_values=simData_liver.variablevalues(:,n-2);
    index_int=ismembertol(SIMTIME,EXPDATA.time{n},1e-4);

    tmpError = tmpError + sum(((simData_values(index_int)- EXPDATA.mean{n}).^2)./(EXPDATA.SD{n}).^2);

end

%% Insulin with only hepatocytes

if n==4% 4: Measured insulin (variable 2)

    simData_values=simData_liver_insulin.variablevalues(:,n-2);
    index_int=ismembertol(SIMTIME,EXPDATA.time{n},1e-4);

    tmpError = tmpError + sum(((simData_values(index_int)- EXPDATA.mean{n}).^2)./(EXPDATA.SD{n}).^2);

end
end

%% Save parameter values that pass the chi-2 test

if (length(tmpError) > 1)
        error = inf;
    else
        if tmpError < CUTOFF
            
            st = '%10.10f ';
            for i = 1 : length(Optpar)
                st = [st '%10.10f '];
            end
            
            st = [st '\n'];
        
            % Write the parameter values and the cost
            fprintf(FID,st,[tmpError Optpar']);
        end
        error = tmpError;
end

end
