% Load parameter indexes

parIndex.i_V_m_hep = ismember(pNames,'V_m_hep');
parIndex.i_V_hep = ismember(pNames,'V_hep');
parIndex.i_V_m_islets = ismember(pNames,'V_m_islets');
parIndex.i_Q = ismember(pNames,'Q');
parIndex.i_EGP_hep = ismember(pNames,'EGP_hep');
parIndex.i_U_ii_hep = ismember(pNames,'U_ii_hep');
parIndex.i_S_i = ismember(pNames,'S_i');
parIndex.i_Sigma = ismember(pNames,'Sigma');
parIndex.i_Alpha = ismember(pNames,'Alpha');
parIndex.i_CL_i_hep = ismember(pNames,'CL_i_hep');
parIndex.i_V_sample_hep = ismember(pNames,'V_sample_hep');
parIndex.i_V_sample_islets = ismember(pNames,'V_sample_islets');

% Disease progession parameters
parIndex.i_Vm = ismember(pNames,'Vm');
parIndex.i_km = ismember(pNames,'km');
parIndex.i_d0 = ismember(pNames,'d0');
parIndex.i_r1 = ismember(pNames,'r1');
parIndex.i_r2 = ismember(pNames,'r2');
parIndex.i_tao_slow = ismember(pNames,'tao_slow');
parIndex.i_G_healthy = ismember(pNames,'G_healthy');
parIndex.i_hep = ismember(pNames,'hep');
parIndex.i_islets = ismember(pNames,'islets');
parIndex.i_k = ismember(pNames,'k');
parIndex.i_kv = ismember(pNames,'kv');
parIndex.i_Vm_f_add = ismember(pNames,'Vm_f_add');
parIndex.i_Km_f_add = ismember(pNames,'Km_f_add');
parIndex.i_Vm_f_islets = ismember(pNames,'Vm_f_islets');
parIndex.i_Km_f_islets = ismember(pNames,'Km_f_islets');

% Media change parameters
parIndex.i_G0 = ismember(pNames,'G0');
parIndex.i_I0 = ismember(pNames,'I0');
parIndex.i_G_GTT = ismember(pNames,'G_GTT');
parIndex.i_I_GTT = ismember(pNames,'I_GTT');
parIndex.i_I_GTT_13 = ismember(pNames,'I_GTT_13');

parIndex.i_delta_G_1 = ismember(pNames,'delta_G_1');
parIndex.i_delta_G_7 = ismember(pNames,'delta_G_7');

parIndex.i_delta_I_1 = ismember(pNames,'delta_I_1');
parIndex.i_delta_I_7 = ismember(pNames,'delta_I_7');

% Media change parameters for each concentration

% 11 mM

parIndex.i_delta_G_1_11mM = ismember(pNames,'delta_G_1_11mM');
parIndex.i_delta_G_7_11mM = ismember(pNames,'delta_G_7_11mM');

parIndex.i_delta_I_1_11mM = ismember(pNames,'delta_I_1_11mM');
parIndex.i_delta_I_7_11mM = ismember(pNames,'delta_I_7_11mM');

% 5.5 mM

parIndex.i_delta_G_7_5p5mM = ismember(pNames,'delta_G_7_5p5mM');
parIndex.i_delta_I_7_5p5mM = ismember(pNames,'delta_I_7_5p5mM');