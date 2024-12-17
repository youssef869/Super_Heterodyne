function BandPassFilt_IF = IF_BPF(signal_num, BW_arr,new_fs)

% RF BPF

fc = 100e3 + signal_num * (50e3);
BW = BW_arr(1,signal_num + 1);

% IF BPF

% Define Bandpass Filter Specifications	
IF_F_pass1 = 25e3 - BW;	
IF_F_pass2 = 25e3 + BW;	
IF_F_stop1 = IF_F_pass1 - 10e3;
IF_F_stop2 = IF_F_pass2 + 10e3;
IF_A_stop1 = 60;	
IF_A_stop2 = 60;		
IF_A_pass = 1;

BandPassSpecObj_IF =  fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
		IF_F_stop1, IF_F_pass1, IF_F_pass2, IF_F_stop2, IF_A_stop1, IF_A_pass, ...
		IF_A_stop2, new_fs);

BandPassFilt_IF = design(BandPassSpecObj_IF, 'butter');

% --- Plot the Filter Characteristics ---
%fvtool(BandPassFilt_IF); % Visualize the filter's frequency and phase response

end