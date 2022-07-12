% Calculates an output energy .mat file from sort_energetics
[times, outputs] = get_output_times();
params = spins_params;

for ii = outputs'
    output_energies = sort_energetics(ii, [0 params.Lx], true);
    offline_diagnos.KE(ii+1) = output_energies.KE_Total;
    offline_diagnos.APE(ii+1) = output_energies.APE_Total;
    offline_diagnos.PE(ii+1) = output_energies.PE_Total;
    offline_diagnos.BPE(ii+1) = output_energies.BPE_Total;
    completion(ii, outputs(end), .1, 'Hillsort');
end

offline_diagnos.Times = times;
save('offline_diagnos.mat', 'offline_diagnos');
