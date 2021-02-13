lib_path='D:\root\LABLOG\FAU\11_CESTpulseq_standard\pulseq-cest-library';

seq_filename=fullfile(lib_path,'/seq-library/APTw_3T_003_2uT_8block_DC95_834ms_braintumor/APTw_3T_003_2uT_8block_DC95_834ms_braintumor.seq')
% read the .seq-file
seq = mr.Sequence;
seq.read(seq_filename);
% get the definitions in the file
offsets_ppm = seq.definitions('offsets_ppm'); % offsets
m0_offset = seq.definitions('M0_offset');     % m0 offset frequency

%% call standard sim
figure(1);
% M_z = Run_pulseq_cest_Simulation(seq_filename,fullfile(lib_path,'/sim-library/GM_3T_001_bmsim.yaml'));
% Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset)
% 
 M_z = Run_pulseq_cest_Simulation(seq_filename,[lib_path '/sim-library/GLIO_3T_001_bmsim.yaml']);
% 
 Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset)

 M_z = Run_pulseq_cest_Simulation(seq_filename,[lib_path '/sim-library/WM_3T_002_bmsim.yaml']);

 Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset)
