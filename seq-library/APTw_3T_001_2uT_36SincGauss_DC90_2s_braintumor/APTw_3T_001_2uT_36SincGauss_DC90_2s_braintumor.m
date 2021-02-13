%% APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor
% Creates a sequence file for an APTw protocol with Sinc-Gaussian pulses, 90% DC and tsat of 2 s
%
% Kai Herz 2020
% kai.herz@tuebingen.mpg.de

% author name for sequence file
author = 'Kai Herz';

%% get id of generation file
if strcmp(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    [~, seqid] = fileparts(which(mfilename));
end

%% sequence definitions
% everything in seq_defs gets written as definition in .seq-file
seq_defs.n_pulses      = 36              ; % number of pulses
seq_defs.tp            = 50e-3           ; % pulse duration [s]
seq_defs.td            = 5e-3            ; % interpulse delay [s]
seq_defs.Trec          = 3.5             ; % recovery time [s]
seq_defs.Trec_M0       = 3.5             ; % recovery time before M0 [s]
seq_defs.M0_offset     = -1560           ; % m0 offset [ppm]
seq_defs.DCsat         = (seq_defs.tp)/(seq_defs.tp+seq_defs.td); % duty cycle
seq_defs.offsets_ppm   = [seq_defs.M0_offset -4:0.25:4]; % offset vector [ppm]
seq_defs.num_meas      = numel(seq_defs.offsets_ppm)   ; % number of repetition
seq_defs.Tsat          = seq_defs.n_pulses*(seq_defs.tp+seq_defs.td) - ...
                         seq_defs.td ;  % saturation time [s]
seq_defs.B0            = 3               ; % B0 [T]
seq_defs.seq_id_string = seqid           ; % unique seq id


%% get info from struct
offsets_ppm = seq_defs.offsets_ppm; % [ppm]
Trec        = seq_defs.Trec;        % recovery time between scans [s]
Trec_M0     = seq_defs.Trec_M0;     % recovery time before m0 scan [s]
tp          = seq_defs.tp;          % sat pulse duration [s]
td          = seq_defs.td;          % delay between pulses [s]
n_pulses    = seq_defs.n_pulses;    % number of sat pulses per measurement. if DC changes use: n_pulses = round(2/(t_p+t_d))
B0          = seq_defs.B0;          % B0 [T]
B1pa        = 1.78;  % mean sat pulse b1 [uT]
spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename

%% scanner limits
% see pulseq doc for more ino
lims = Get_scanner_limits();

%% create scanner events
% satpulse
gyroRatio_hz  = 42.5764;                  % for H [Hz/uT]
gyroRatio_rad = gyroRatio_hz*2*pi;        % [rad/uT]
fa_sat        = B1pa*gyroRatio_rad*tp; % flip angle of sat pulse
% create pulseq saturation pulse object

%satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', t_p,'system',lims,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss
satPulse      = mr.makeSincPulse(fa_sat, 'Duration', tp, 'system', lims,'timeBwProduct', 2,'apodization', 0.15); % philips-like sinc

[B1cwpe,B1cwae,B1cwae_pure,alpha]= calc_power_equivalents(satPulse,tp,td,1,gyroRatio_hz);
seq_defs.B1cwpe = B1cwpe;

% spoilers
spoilRiseTime = 1e-3;
spoilDuration = 4500e-6+ spoilRiseTime; % [s]
% create pulseq gradient object
[gxSpoil, gySpoil, gzSpoil] = Create_spoiler_gradients(lims, spoilDuration, spoilRiseTime);

% pseudo adc, not played out
pseudoADC = mr.makeAdc(1,'Duration', 1e-3);

%% loop through zspec offsets
offsets_Hz = offsets_ppm*gyroRatio_hz*B0;

% init sequence
seq = mr.Sequence();
% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    if currentOffset == seq_defs.M0_offset*gyroRatio_hz*B0
        if Trec_M0 > 0
            seq.addBlock(mr.makeDelay(Trec_M0));
        end
    else
        if Trec > 0
            seq.addBlock(mr.makeDelay(Trec)); % recovery time
        end
    end
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
    accumPhase=0;
    for np = 1:n_pulses
        satPulse.phaseOffset = mod(accumPhase,2*pi); % set accumulated pahse from previous rf pulse
        seq.addBlock(satPulse) % add sat pulse
        % calc phase for next rf pulse
        accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
        if np < n_pulses % delay between pulses
            seq.addBlock(mr.makeDelay(td)); % add delay
        end
    end
    if spoiling % spoiling before readout
        seq.addBlock(gxSpoil,gySpoil,gzSpoil);
    end
    seq.addBlock(pseudoADC); % readout trigger event
end



%% write definitions
def_fields = fieldnames(seq_defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

%% plot
save_seq_plot(seq_filename);

%% call standard sim
M_z = Run_pulseq_cest_Simulation(seq_filename,'../../sim-library/GM_3T_001_bmsim.yaml');

%% normalization
% The normalization depends on the M0 offset an their positions in the .seq-file
% In this example, the M0 offset is set as a definition in the APTw_3T_example.seq
% It can be accessed with the correct key and is used as the normalization
% value

% read the .seq-file
seq = mr.Sequence;
seq.read(seq_filename);

% get the definitions in the file
offsets_ppm = seq.definitions('offsets_ppm'); % offsets
m0_offset = seq.definitions('M0_offset');     % m0 offset frequency

% we search for all appearances of the m0 offset and normalize by it, or by
% the mean value if all m0 offsets if multiple occur
M0_idx = [];
if any(offsets_ppm == m0_offset)
    M0_idx = [M0_idx find(offsets_ppm == m0_offset)];
    M0 = mean(M_z(M0_idx));
    M_z(M0_idx) = [];
    offsets_ppm(M0_idx) = [];
    Z = M_z./M0; % calculate the Z-value
end

%% plot Z-spectrum and MTRasym
% sort normalized values 
[ppm_sort, idx] = sort(offsets_ppm);
Z = Z(idx);
% MTRasym -> for this example we have the same distribution of offsets on
% both sides of the water, so we can calculate the asym easily.
% For different .seq-files, different steps may be neccessary
MTRasym=Z(end:-1:1)-Z;
MTRasym(1:ceil(end/2)) = NaN;

% finally, plot the Z-spectrum and the MTRasym curve
figure; hold on;
yyaxis left;
plot(ppm_sort, Z,'Displayname','Z-spectrum'); set(gca,'xdir','reverse');
xlabel('\Delta\omega [ppm]'); legend show;
yyaxis right;
plot(ppm_sort,MTRasym,'Displayname','MTR_{asym}');
axis([ppm_sort(1) ppm_sort(end) -(2*max(abs(MTRasym))) (2*max(abs(MTRasym)))])





