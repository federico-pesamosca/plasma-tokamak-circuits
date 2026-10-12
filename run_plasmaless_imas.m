% run the plasmaless coil-vessel model through its IMAS interface
% (same inputs as run_plasmaless_model.m): the model is built from
% machine-description IDSs, and every time step takes a one-slice pf_active
% with the coil voltages at t and returns pf_active/pf_passive currents at t+dt

clearvars -except do_plot
% plot coil currents at the end (set do_plot = 1 before calling to enable).
% Default 0: figures hang with MATLAB R2025b on the SDCC login nodes; use
% MATLAB R2024b (module load MATLAB/2024b-r5) to plot with plot_plasmaless_model.
if ~exist('do_plot','var'), do_plot = 0; end

% IMAS interface functions in imas_model/ (next to this script)
repo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(repo_root,'imas_model'));

%% simulation parameters
tstart = 0;
tend = 28;
dt = 0.001;
write_output = 0;   % 1: write the output slices to out_uri
out_uri = 'imas:hdf5?path=/scratch/users/schneim/plasmaless/run_out';
results_file = fullfile(repo_root,'plasmaless_imas_results.mat');  % saved time, Ia, Ie, names

%% load machine description
% DD4 entry shipped in data/md_dd4, made with tools/convert_md_dd3_to_dd4.py
% (IMAS-MATLAB DD-4.x module)
md_uri = ['imas:hdf5?path=' fullfile(repo_root,'data','md_dd4')];
fprintf('Loading machine description from %s\n', md_uri);
[em_coupling,pf_active_md,pf_passive_md] = load_plasmaless_md(md_uri);

%% build model and initial state (cold start)
[model, x] = configure_plasmaless_model(em_coupling, pf_active_md, pf_passive_md, dt);
fprintf('Model built: %d coils, %d passive loops, dt = %g s\n', model.n_coils, model.n_loops, dt);
% warm start from a one-slice pf_active/pf_passive pair:
% [model, x] = configure_plasmaless_model(em_coupling, pf_active_md, pf_passive_md, dt, ...
%     pf_active_slice, pf_passive_slice);

%% prepare inputs
TS = timeseries_plasmaless_model(tstart,tend,dt);

%% run model (one IMAS step per time step, see imas_model/simulator_plasmaless_model_imas.m)
if write_output, uri_arg = out_uri; else, uri_arg = []; end
[time, Ia, Ie] = simulator_plasmaless_model_imas(TS, model, x, uri_arg);

%% save results and plot (plot_plasmaless_model.m)
coil_names = model.coil_names;
loop_names = model.loop_names;
save(results_file,'time','Ia','Ie','coil_names','loop_names','-v7');
fprintf('Results saved to %s\n', results_file);

% same output structure as the Simulink/MATLAB simulators of run_plasmaless_model.m:
% out.simout is a timeseries of [Ia Ie] (nt x (n_coils+n_loops)); TS, tstart, tend as above
out.tout = time;
out.simout = timeseries([Ia Ie], time);
if do_plot, plot_plasmaless_model; end
