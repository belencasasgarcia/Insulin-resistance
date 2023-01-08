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

%% Simulation for the co-culture, hyperglycemia

Optparam_11mM=Optpar;
Optparam_11mM(parIndex.i_G0)=11;
Optparam_11mM(parIndex.i_I_GTT)=0;
Optparam_11mM(parIndex.i_delta_G_1)=Optparam_11mM(parIndex.i_delta_G_1_11mM);
Optparam_11mM(parIndex.i_delta_G_13)=Optparam_11mM(parIndex.i_delta_G_13_11mM);
Optparam_11mM(parIndex.i_delta_I_1)=Optparam_11mM(parIndex.i_delta_I_1_11mM);
Optparam_11mM(parIndex.i_delta_I_13)=Optparam_11mM(parIndex.i_delta_I_13_11mM);

icOrig_11mM=[(11+Optparam_11mM(parIndex.i_delta_G_1))*Optparam_11mM(parIndex.i_V_m_hep) ...
    (11+Optparam_11mM(parIndex.i_delta_G_1))*Optparam_11mM(parIndex.i_V_m_islets)...
   (0+Optparam_11mM(parIndex.i_delta_I_1))*Optparam_11mM(parIndex.i_V_m_hep) (0+Optparam_11mM(parIndex.i_delta_I_1))*Optpar(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088 0 0];

try
    simData_islets_liver = feval(char(modelName), SIMTIME, icOrig_11mM, Optparam_11mM, COSTOPTIONS);
    
catch simFail
    
    disp('Simulation crashed');
    error = inf;
    disp(param);
    disp(simFail)
    return;
end

%% Simulation for 5.5 mM dose

Optparam_5p5mM=Optpar;
Optparam_5p5mM(parIndex.i_G0)=5.5;
Optparam_5p5mM(parIndex.i_I_GTT)=0;
Optparam_5p5mM(parIndex.i_delta_G_1)=Optparam_5p5mM(parIndex.i_delta_G_1_5p5mM);
Optparam_5p5mM(parIndex.i_delta_G_13)=Optparam_5p5mM(parIndex.i_delta_G_13_5p5mM);
Optparam_5p5mM(parIndex.i_delta_I_1)=Optparam_5p5mM(parIndex.i_delta_I_1_5p5mM);
Optparam_5p5mM(parIndex.i_delta_I_13)=Optparam_5p5mM(parIndex.i_delta_I_13_5p5mM);

icOrig_5p5=[(5.5+Optparam_5p5mM(parIndex.i_delta_G_1))*Optparam_5p5mM(parIndex.i_V_m_hep) ...
    (11+Optparam_5p5mM(parIndex.i_delta_G_1))*Optparam_5p5mM(parIndex.i_V_m_islets)...
   (0+Optparam_5p5mM(parIndex.i_delta_I_1))*Optparam_5p5mM(parIndex.i_V_m_hep) (0+Optparam_5p5mM(parIndex.i_delta_I_1))*Optpar(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088 0 0];

try
    simData_islets_liver_5p5 = feval(char(modelName), SIMTIME, icOrig_5p5, Optparam_5p5mM, COSTOPTIONS);
    
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
% 5: Measured glucose with liver+islets, 5.5 mM   (variable 1)
% 6: Measured insulin with liver+islets, 5.5 mM   (variable 2)


if n==1 || n==2
    
    % Values from the co-culture simulation
    simData_values=simData_islets_liver.variablevalues(:,n);
    index_int=ismembertol(SIMTIME,EXPDATA.time{n});
    
    tmpError = tmpError + sum(((simData_values(index_int)- EXPDATA.mean{n}).^2)./(EXPDATA.SD{n}).^2);
    
% figure()
% plot(SIMTIME(index_int),simData_values(index_int))
% hold on
% plot(EXPDATA.time{n},EXPDATA.mean{n})

end

%% Liver and islets, 5.5 mM 

if n==5 || n==6
    
    % Values from the co-culture simulation
    simData_values=simData_islets_liver_5p5.variablevalues(:,n-4);
    index_int=ismembertol(SIMTIME,EXPDATA.time{n});

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
