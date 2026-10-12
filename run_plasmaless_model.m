% run simulator_plasmaless.slx (simulink) 
% or simulator_plasmaless_model_matlab.m (matlab)
% featuring dynamical plasmaless coil-vessel circuit model
% and plot the results

clearvars -except do_plot
% plot at the end (set do_plot = 1 before calling to enable). Default 0:
% plotting in R2025b hangs on SDCC login nodes; use R2024b (plot_plasmaless_model)
if ~exist('do_plot','var'), do_plot = 0; end

%% flag for executing matlab or simulink mode
mode = 'simulink'; % 'matlab' or 'simulink'

%% simulation parameters
tstart = 0;
tend = 28;
dt = 0.001;

%% load model (DD4 machine description in data/md_dd4)
[em_coupling,pf_active,pf_passive] = load_plasmaless_md();
[model,x0] = configure_plasmaless_model(em_coupling,pf_active,pf_passive,dt);
% matrices used by the simulink model and by the matlab simulator (base workspace)

%% prepare timeseries
TS = timeseries_plasmaless_model(tstart,tend,dt);

%% run model
switch mode
    case 'simulink'
        A = model.A; B = model.B; C = model.C; D = model.D;
        out = sim('simulator_plasmaless_model_simulink.slx');
    case 'matlab'
        out = simulator_plasmaless_model_matlab(TS,model,x0);
end

%% save results (same variables and units as run_plasmaless_imas.m)
time = out.simout.Time;
y = squeeze(out.simout.Data);
if size(y,1) ~= numel(time), y = y.'; end
n_coils = model.n_coils;
Ia = y(:,1:n_coils);
Ie = y(:,n_coils+1:end);
coil_names = model.coil_names;
loop_names = model.loop_names;
results_file = fullfile(fileparts(mfilename('fullpath')),'plasmaless_matlab_results.mat');
save(results_file,'time','Ia','Ie','coil_names','loop_names','-v7');
fprintf('Results saved to %s\n', results_file);

%% plot model
if do_plot, plot_plasmaless_model; end