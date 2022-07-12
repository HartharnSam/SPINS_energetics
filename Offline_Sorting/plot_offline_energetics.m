%PLOT_OFFLINE_ENERGETICS
%
% Recalculates energetics from Offline SPINS diagnostics, by:
% Getting the energy partitioned into KE, PE, BPE, APE
% Then calculates the derivatives dKE_dt, dPE_dt, dBPE_dt, dAPE_dt
% Then calculates the pathways of the changes (phi_i, phi_m, dissipation,
% etc)
%
% See also plot_diagnos, calc_energetics

clc; clearvars; close all;
% Optional switch to compare to SPINS's online diagnostics / plot_diagnos
% output
compareSPINS = true;

if compareSPINS
    load all_diagnos
    % Set up temporal grid
    time = all_diagnos.diagnos.Time;
end
%% Load in offline diagnostics from hill_sort
m_path = mpath;
addpath([m_path, '../01_sort_hill/'])
times = first_output():last_output();
recalculate = true;

if recalculate
    run_sort_energetics;
else
    load('offline_diagnos');
    disp('Loaded (not recalculated) offline diagnostics')
end

hillsort_dt = times(2)-times(1);
max_time = max(times);

%% First look at energy components
fig = figure(1);
n_plots_1 = 5;
% Potential Energy
subaxis(n_plots_1, 1, 1)
plot(times, offline_diagnos.PE, 'bx')
if compareSPINS
hold on
PE = all_diagnos.diagnos.PE_tot;
plot(time, PE, 'k-');

end
xlim([0 max_time]);
ylabel('PE Total')

% BPE
subaxis(n_plots_1, 1, 2);
plot(times, (offline_diagnos.BPE-offline_diagnos.BPE(1))/2.3, 'bx'); % TODO: Resolve this dodgy af fix
if compareSPINS
    hold on
    BPE = all_diagnos.diagnos.BPE_tot;
    plot(time, BPE-BPE(1), 'k-');
end
xlim([0 max_time]);
ylabel('BPE tot');

% APE
subaxis(n_plots_1, 1, 3);
APE = PE - BPE;
plot(time, APE, '-', 'Color', [.3 .3 .3]);
hold on
plot(times, offline_diagnos.APE-11.9, 'bx'); % TODO: Resolve this dodgy af fix
xlim([0 max_time]);
ylabel('APE tot');

% KE
subaxis(n_plots_1, 1, 4);
plot(times, offline_diagnos.KE, 'bx');
if compareSPINS
    KE = all_diagnos.diagnos.KE_x + all_diagnos.diagnos.KE_z;
    plot(time, KE, 'k-');
    hold on
end
xlim([0 max_time]);
ylabel('KE');
legend('SPINS', 'MATLAB');

% Total Energy
subaxis(n_plots_1, 1, 5)
%E_tot = all_diagnos.EnergyBudget.E_tot;
TotE = offline_diagnos.KE+offline_diagnos.PE;
plot(times, TotE, 'bx');
xlim([0 max_time]);
ylabel('E Total');

figure_print_format(gcf)
fig = gcf;
fig.Position([3 4]) = [768 432.8000];


%% Now look at Energy Rates
fig2 = figure(2);
n_plots_2 = 4;
% Get the finite differencing matrix
Dmat = FiniteDiff(times, 1, 2, true);

%dPE_dt
subaxis(n_plots_2, 1, 1);
dPE_dt = Dmat*offline_diagnos.PE';
plot(times, dPE_dt, 'bx');
xlim([0 max_time]);
ylabel('dPE/dt');

%dBPE_dt
subaxis(n_plots_2, 1, 2);
dBPE_dt = Dmat*offline_diagnos.BPE';
plot(times, dBPE_dt, 'bx');
xlim([0 max_time]);
ylabel('dBPE/dt');

%dAPE_dt
subaxis(n_plots_2, 1, 3);
dAPE_dt = Dmat*offline_diagnos.APE';
plot(times, dAPE_dt, 'bx');
if compareSPINS
    hold on
    plot(all_diagnos.EnergyRates.Time, all_diagnos.EnergyRates.APE_rate, 'b--')
end
xlim([0 max_time])
ylabel('dAPE/dt')

%dKE_dt
subaxis(n_plots_2, 1, 4);
dKE_dt = Dmat*offline_diagnos.KE';
plot(times, dKE_dt, 'bx')
if compareSPINS
    hold on
    plot(all_diagnos.EnergyRates.Time, all_diagnos.EnergyRates.KE_rate, 'b--')
end
xlim([0 max_time])
ylabel('dKE/dt');

figure_print_format(fig2)
fig2.Position([3 4]) = [768 432.8000];

%% Conversions

