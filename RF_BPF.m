function BandPassFilt_RF = RF_BPF(signal_num, BW_arr,new_fs)

% RF BPF

fc = 100e3 + signal_num * (50e3);
BW = BW_arr(1,signal_num + 1);

% Define Bandpass Filter Specifications
F_pass1 = (fc - BW);	% Edge of the passband 
F_pass2 = (fc + BW);	% Closing edge of the passband
F_stop1 = F_pass1 - 10e3;	% Edge of the stopband
F_stop2 = F_pass2 + 10e3;	% Edge of the second stopband
A_stop1 = 60;		% Attenuation in the first stopband
A_stop2 = 60;		% Attenuation in the second stopband
A_pass = 1;		% Amount of ripple allowed in the passband

BandPassSpecObj_RF =  fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
		F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, ...
		A_stop2, new_fs);

BandPassFilt_RF = design(BandPassSpecObj_RF, 'butter');

% --- Plot the Filter Characteristics ---
fvtool(BandPassFilt_RF); % Visualize the filter's frequency and phase response

end