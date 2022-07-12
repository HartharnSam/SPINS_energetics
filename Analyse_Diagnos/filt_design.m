%FILT_DESIGN - This is just to explain the choice of BPE filter, compared
%to the one used in plot_diagnos. Three filters are explored, mlptdenoise,
%a spectral lowpass filter, and the MATLAB smooth filter (running mean). 
%
% The filter is needed to smooth out tiny variations in BPE, calculated at
% each timestep using the sorting algorithm. To calculate mixing, we need
% dBPE/dt, which then becomes really noisy without smoothing of the BPE
% signal. 

clc; clearvars; close all;

load all_diagnos

BPE = all_diagnos.diagnos.BPE_tot;
Time = all_diagnos.diagnos.Time;

dt = 0.01;
time_rate = Time(1):dt:Time(end);
BPE_const_t = interp1(Time, BPE, time_rate, 'pchip')';

%% Option 1 : mlptdenoise
% This is the filter used by plot_diagnos - I assume there was a good
% reason for this, however it needs 1) the wavelet toolbox, 2) it produces
% an NxN matrix during the calculation which for most simulations gets too
% big for MATLAB to handle 

%tic;
%BPE_tot = mlptdenoise(BPE, Time, 5, 'DualMoments', 3);
%toc;
 
% Some sample computation times
close all; 
n = [2527 5054 6318 10109 12637 25274 50549];
t = [0.7068 2.793921 4.240124 9.249551 14.156055 52.63 0];
plot(n, t);
hold on
plot(n, (n/n(1)).^2);
legend('Calculation time', 'n^{2}');
% We can see the computation time increases to almost n^2, so it isn't
% appropriate to use for long simulations with lots of timesteps. 

%% Option 2 : Lowpass filter
% Possibly the most intuative filter to use here, as we are trying to
% remove high frequency noise. Cut off the abs(k)>k_cutoff wavenumbers. 
% It is fast, it needs a regularly spaced grid

k_cutoff = 30;
fft_y0 = fftshift(fft(BPE_const_t));
fft_y0([1:k_cutoff+1 end-k_cutoff:end]) = 0;
data = real(ifft(ifftshift(fft_y0)));

figure(2)
plot(time_rate, data, 'k-');
% BUT - it does some odd things at the start and end of the timeseries
% making it unmanageable. 
hold on
plot(Time, BPE, 'b-');

%% Option 3 : Smooth filter (running mean)
% Maybe the least justifiable method (its hardest to physically justify our 
% cutoff value), and we can lose info from all scales. 
% But it's easy to use, implement, and super fast to run

smooth_window = 20;
BPE_tot = smooth(BPE_const_t, smooth_window);
plot(time_rate, BPE_tot, 'r-');
legend('raw', 'fft filt', 'smooth filter')

pause;
xlim([5 10])
%Crucially, it also works