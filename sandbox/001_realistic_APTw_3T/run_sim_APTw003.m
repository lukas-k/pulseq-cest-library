lib_path='D:\root\LABLOG\FAU\11_CESTpulseq_standard\pulseq-cest-library';

seq_filename=fullfile(lib_path,'/seq-library/APTw_3T_003_2uT_8block_DC95_834ms_braintumor/APTw_3T_003_2uT_8block_DC95_834ms_braintumor.seq')
% read the .seq-file
seq = mr.Sequence;
seq.read(seq_filename);
% get the definitions in the file
offsets_ppm = seq.definitions('offsets_ppm'); % offsets
m0_offset = seq.definitions('M0_offset');     % m0 offset frequency

%% call standard sim  for WM and GLIO ( with shifted Lorentzian MT)
figure('Name','(Super) Lorentz');
% M_z = Run_pulseq_cest_Simulation(seq_filename,fullfile(lib_path,'/sim-library/GM_3T_001_bmsim.yaml'));
% Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset)
% 
sandbox_path='D:\root\LABLOG\FAU\11_CESTpulseq_standard\pulseq-cest-library\sandbox\001_realistic_APTw_3T';
 M_z = Run_pulseq_cest_Simulation(seq_filename,[sandbox_path '/WM_3T_001_bmsim.yaml']);
 [Z,w]= Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset);

 M_z = Run_pulseq_cest_Simulation(seq_filename,[sandbox_path '/WM_3T_003_bmsim.yaml']);
 [Z1,w]=Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset);
 
%% call standard sim  for WM and GLIO ( with shifted Lorentzian MT)
figure('Name','MTC asym');
% M_z = Run_pulseq_cest_Simulation(seq_filename,fullfile(lib_path,'/sim-library/GM_3T_001_bmsim.yaml'));
% Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset)
% 
sandbox_path='D:\root\LABLOG\FAU\11_CESTpulseq_standard\pulseq-cest-library\sandbox\001_realistic_APTw_3T';
 M_z = Run_pulseq_cest_Simulation(seq_filename,[sandbox_path '/WM_3T_003_bmsim.yaml']);
 [Z,w]= Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset);

 M_z = Run_pulseq_cest_Simulation(seq_filename,[sandbox_path '/GLIO_3T_003_bmsim.yaml']);
 [Z1,w]=Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset);
 
 %% call standard sim  for WM and GLIO ( with shifted Lorentzian MT)
figure('Name','WM vs GLIO (same CEST pools)');
% M_z = Run_pulseq_cest_Simulation(seq_filename,fullfile(lib_path,'/sim-library/GM_3T_001_bmsim.yaml'));
% Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset)
% 
sandbox_path='D:\root\LABLOG\FAU\11_CESTpulseq_standard\pulseq-cest-library\sandbox\001_realistic_APTw_3T';
 M_z = Run_pulseq_cest_Simulation(seq_filename,[sandbox_path '/WM_3T_002_bmsim.yaml']);
 [Z,w]= Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset);

 M_z = Run_pulseq_cest_Simulation(seq_filename,[sandbox_path '/GLIO_3T_001_bmsim.yaml']);
 [Z1,w]=Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset);
 
  %% call standard sim  for WM and GLIO ( with shifted Lorentzian MT)
figure('Name','WM vs GLIO (low GLIO_NOE)');
% M_z = Run_pulseq_cest_Simulation(seq_filename,fullfile(lib_path,'/sim-library/GM_3T_001_bmsim.yaml'));
% Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset)
% 
sandbox_path='D:\root\LABLOG\FAU\11_CESTpulseq_standard\pulseq-cest-library\sandbox\001_realistic_APTw_3T';
 M_z = Run_pulseq_cest_Simulation(seq_filename,[sandbox_path '/WM_3T_002_bmsim.yaml']);
 [Z,w]= Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset);

 M_z = Run_pulseq_cest_Simulation(seq_filename,[sandbox_path '/GLIO_3T_002_bmsim.yaml']);
 [Z1,w]=Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset);
 
  %% EMR?
figure('Name','EMR in glioma');
  M_z = Run_pulseq_cest_Simulation(seq_filename,[sandbox_path '/GLIO_3T_002_bmsim.yaml']);
 [Z,w]= Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset);

 M_z = Run_pulseq_cest_Simulation(seq_filename,[sandbox_path '/GLIO_3T_003_bmsim.yaml']);
 [Z1,w]=Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset);
  
 figure(3), plot(w,(Z1-Z)); set(gca,'xdir','reverse');
 
 
 
