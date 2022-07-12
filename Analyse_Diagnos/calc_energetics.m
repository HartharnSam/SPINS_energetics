% Recalculates energetics from online SPINS diagnostics, by:
% Getting the energy partitioned into KE, PE, BPE, APE
% Then calculates the derivatives dKE_dt, dPE_dt, dBPE_dt, dAPE_dt
% Then calculates the pathways of the changes (phi_i, phi_m, dissipation,
% etc)
% For each, it compares to the data produced by plot_diagnos;
% Dev TODO: 
%   - Totals for rates (KE2APE, KE2Int, APE2BPE, Int2BPE)
%   - Do we want to output scales in here? Mixing efficiencies?
%   - end dependence on the plot_diagnos script?

clc; clearvars; close all;
compareSPINS = true;

% Load in spins diagnostics
load all_diagnos
time = all_diagnos.diagnos.Time;
max_time = max(time);

% Make a regular time grid, with a pre-determined dt (don't set it too
% small like it is in plot_diagnos)
dt = 0.01;
time_rate = time(1):dt:time(end);
%Alternatively the definition in the original plot_diagnos:
% time_rate = linspace(time(1), time(end), round(length(time)/2)); % Linear spaced simulation times
% dt = time_rate(2)-time_rate(1);

%% First look at energy components
fig1 = figure(1);
n_plots_1 = 7;
% Potential Energy
subaxis(n_plots_1, 1, 1)
PE = all_diagnos.diagnos.PE_tot;
% Interpolate this onto a constant time grid like in plot_diagnos (L319)
PE_const_t = interp1(time, PE, time_rate, 'pchip')';
plot(time_rate, PE_const_t, 'k-');
if compareSPINS
    hold on
    plot(time, PE, 'b--');
end
xlim([0 max_time]);
ylabel('PE Total')

% BPE
subaxis(n_plots_1, 1, 2);
BPE = all_diagnos.diagnos.BPE_tot;
BPE_const_t = interp1(time, BPE, time_rate, 'pchip')';
% BPE needs a smoothing or filter. See different filters decisions in the
% filt_design.m script
%k_cutoff = 40;
%fft_y0 = fftshift(fft(BPE_const_t));
%fft_y0([1:k_cutoff+1 end-k_cutoff:end]) = 0;
%BPE_const_t = real(ifft(ifftshift(fft_y0)));
BPE_const_t = smooth(BPE_const_t, 10);
plot(time_rate, BPE_const_t, 'k-');
if compareSPINS
    hold on
    plot(time, BPE, 'b--');
end
xlim([0 max_time]);
ylabel('BPE tot');

% APE
subaxis(n_plots_1, 1, 3);
APE = PE - BPE;
APE_const_t = PE_const_t - BPE_const_t;
plot(time_rate, APE_const_t, 'k-');
if compareSPINS
    hold on
    plot(time, APE, 'b--');
end
xlim([0 max_time]);
ylabel('APE tot');

% KE
subaxis(n_plots_1, 1, 4);
KE = all_diagnos.diagnos.KE_x + all_diagnos.diagnos.KE_z;
KE_const_t = interp1(time, KE, time_rate, 'pchip')';
plot(time_rate, KE_const_t, 'k-');
if compareSPINS
    hold on
    plot(time, KE, 'b--');
end
xlim([0 max_time]);
ylabel('KE');
legend('Calculated', 'SPINS');

% Total Energy
subaxis(n_plots_1, 1, 5)
E_tot = all_diagnos.EnergyBudget.E_tot;
TotE = KE_const_t+PE_const_t;
plot(time_rate, TotE, 'k-');
if compareSPINS
    hold on
    plot(all_diagnos.EnergyBudget.Time, E_tot, 'b--');
end
xlim([0 max_time]);
ylabel('E Total');

figure_print_format(fig1)
monitor_loc = groot().MonitorPositions;
num_monitors = size(monitor_loc, 1);
if num_monitors == 2
    fig1.Position = [1924 50.6 768 729];
else
    fig1.Position = [5.8000 217 536.8000 558.4000];
end

%% Now look at Energy Rates
fig2 = figure(2);
n_plots_2 = 6;
% Get the finite differencing matrix
Dmat = FiniteDiff(time_rate, 1, 2, true);

%dPE_dt
subaxis(n_plots_2, 1, 1);
dPE_dt = Dmat*PE_const_t;
plot(time_rate(10:end-10), dPE_dt(10:end-10), 'k-');
xlim([0 max_time]);
ylabel('dPE/dt');

% dBPE_dt
subaxis(n_plots_2, 1, 2);
dBPE_dt = Dmat*BPE_const_t;
% We could re-filter dBPE_dt?
%k_cutoff = length(time_rate)/3; % This is SOOOO aggressive and it still doesn't really denoise
%fft_y0 = fftshift(fft(dBPE_dt));
%fft_y0([1:k_cutoff+1 end-k_cutoff:end]) = 0;
%dBPE_dt = real(ifft(ifftshift(fft_y0)));
%dBPE_dt = smooth(dBPE_dt, 20);

plot(time_rate(10:end-10), dBPE_dt(10:end-10), 'k-');
if compareSPINS
    % You probably really don't want this mess ;-)
    %hold on
    %plot(all_diagnos.EnergyRates.Time, all_diagnos.EnergyRates.BPE_rate, 'b--')
end
xlim([0 max_time]);
ylabel('dBPE/dt');

% dAPE_dt
subaxis(n_plots_2, 1, 3);
dAPE_dt = Dmat*APE_const_t;
% % We could re-filter dAPE_dt? %FIXME
% k_cutoff = length(time_rate)/3; % This is SOOOO aggressive and it still doesn't really denoise
% fft_y0 = fftshift(fft(dAPE_dt));
% fft_y0([1:k_cutoff+1 end-k_cutoff:end]) = 0;
% dAPE_dt = real(ifft(ifftshift(fft_y0)));
%dAPE_dt = smooth(dAPE_dt, 20);

plot(time_rate(10:end-10), dAPE_dt(10:end-10), 'k-');
if compareSPINS
    hold on
    plot(all_diagnos.EnergyRates.Time, all_diagnos.EnergyRates.APE_rate, 'b--')
end
xlim([0 max_time]);
ylabel('dAPE/dt');

% dKE_dt
subaxis(n_plots_2, 1, 4);
dKE_dt = Dmat*KE_const_t;
plot(time_rate(10:end), dKE_dt(10:end), 'k-');
if compareSPINS
    hold on
    plot(all_diagnos.EnergyRates.Time, all_diagnos.EnergyRates.KE_rate, 'b--')
end
xlim([0 max_time]);
ylabel('dKE/dt');

% dInt_dt
subaxis(n_plots_2, 1, 5);
Diss = all_diagnos.diagnos.Diss_tot;
Diss_const_t = interp1(time, Diss, time_rate, 'pchip')';

Int2BPE = all_diagnos.diagnos.BPE_from_int;
phi_i_const_t = interp1(time, Int2BPE, time_rate, 'pchip')';

dInt_dt = Diss_const_t - phi_i_const_t;

plot(time_rate(10:end), dInt_dt(10:end), 'k-');
xlim([0 max_time]);
ylabel('dInt/dt');

figure_print_format(fig2)
if num_monitors == 2
    fig2.Position = [2.5746e+03 50.6000 599.2000 729];
else
    fig2.Position = [464.2000 342.6000 592 432.8000];
end

%% Now look at conversions
fig3 = figure(3);
n_plots_3 = 4;
subaxis(n_plots_3, 1, 1);
plot(time_rate(10:end-10), Diss_const_t(10:end-10), 'k-');
if compareSPINS
    hold on
    plot(time(10:end), Diss(10:end), 'b--');
end
xlim([0 max_time]);
ylabel('$\epsilon$', 'interpreter', 'latex');

subaxis(n_plots_3, 1, 2);
KE2APE = all_diagnos.EnergyRates.KE2APE_rate;
phi_z = -dKE_dt - Diss_const_t; %KE2APE
plot(time_rate(10:end-10), phi_z(10:end-10), 'k-');
% if compareSPINS % These are plot_diagnos outputs rather than from SPINS
%     hold on
%     plot(all_diagnos.EnergyRates.Time(10:end), KE2APE(10:end), 'b--');
% end
xlim([0 max_time]);
ylabel('$\phi_z$', 'interpreter', 'latex');

subaxis(n_plots_3, 1, 3);
APE2BPE = all_diagnos.EnergyRates.APE2BPE_rate;
phi_m = phi_z - dAPE_dt; %APE2BPE

phi_m = dBPE_dt - phi_i_const_t;

plot(time_rate(10:end-10), phi_m(10:end-10), 'k-');
% if compareSPINS % plot_diagnos outputs, rather than from spins
%plot(all_diagnos.EnergyRates.Time(10:end), APE2BPE(10:end), 'b--');
%hold on
%end
xlim([0 max_time]);
ylabel('$\phi_m$', 'interpreter', 'latex');

subaxis(n_plots_3, 1, 4)
phi_i = phi_i_const_t;%dBPE_dt - phi_m;
plot(time_rate(10:end-10), phi_i(10:end-10)*1000, 'k-');
if compareSPINS
    hold on
    plot(time(10:end), Int2BPE(10:end)*1000, 'b--');
end
xlim([0 max_time]);
ylabel('$\phi_i$', 'interpreter', 'latex');

figure_print_format(fig3);
if num_monitors == 2
    fig3.Position = [3.1746e+03 349 644 432.8000];
else
    fig3.Position = [464.2000 342.6000 592 432.8000];
end

%% Finally Calculate some things we weren't able to before:

TotE_loss = TotE(1) - TotE;        
% NOTE: Actual Internal Energy isn't ever known, but SPINS does calculate how it
% changes, from this we can calculate a "Internal Energy Anomaly"
figure(1);
subaxis(n_plots_1, 1, 6);
Int_const_t = cumtrapz(time_rate, dInt_dt); 
plot(time_rate, Int_const_t, 'k-');
xlim([0 max_time]);
ylabel('Int. E');

subaxis(n_plots_1, 1, 7);
Num_const_t  = TotE_loss - dInt_dt; % There should also be a forcing term here if forcing is on
plot(time_rate(2:end), Num_const_t(2:end), 'k-');
xlim([0 max_time]);
ylabel('Num. Lost');

figure_print_format(gcf);

figure(2);
subaxis(n_plots_2, 1, 6);
dNumE_dt = Dmat*Num_const_t;
plot(time_rate(10:end-10), dNumE_dt(10:end-10), 'k-');
xlim([0 max_time]);
ylabel('dNum/dt');

figure_print_format(gcf);

%% Put this all into a usable structure
% How much energy has gone into the numerics (tot)?

EnergyBudget = struct('Time', time_rate', 'E_tot', TotE, 'KE_tot', KE_const_t, ...
    'APE_tot', APE_const_t, 'BPE_tot', BPE_const_t, 'PE_tot', PE_const_t);
EnergyRates = struct('Time', time_rate', 'dKE_dt', dKE_dt, 'dAPE_dt', dAPE_dt,...
    'dBPE_dt', dBPE_dt, 'dPE_dt', dPE_dt,...
    'Diss', Diss_const_t, 'APE2BPE', phi_m, 'Int2BPE', phi_i, 'KE2APE', phi_z);

save('energy_diagnos', 'EnergyBudget', 'EnergyRates');
