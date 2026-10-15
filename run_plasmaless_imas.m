% run the plasmaless coil-vessel model through its IMAS interface
% (same inputs as run_plasmaless_model.m): the model is built from
% machine-description IDSs (or an external linear model, see
% linear_model_dir, or the vertical-instability model, see vertical_uri),
% and every time step takes a one-slice pf_active
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
% '' = circuit model built from the machine description; the path of a directory
% holding linear_model.json + matrices = external linear model (deviations from
% its linearisation point), e.g. fullfile(repo_root,'linearModelSimple_105073_t84.417598').
% Default: the environment variable PLASMALESS_LINEAR_MODEL_DIR when set, else ''.
linear_model_dir = getenv('PLASMALESS_LINEAR_MODEL_DIR');
% Vertical-instability model (vertical_instability/, rigid vertical displacement
% around an equilibrium, states = deviations from its currents): set
% vertical_uri to an IMAS URI or entry directory with an equilibrium holding
% profiles_2d. Defaults: the environment variables PLASMALESS_VERTICAL_EQUILIBRIUM_URI
% (else '' = not used), PLASMALESS_VERTICAL_EQUILIBRIUM_TIME (else the extractor's
% default time) and PLASMALESS_VERTICAL_PASSIVE_CURRENTS ('equilibrium' or 'zero',
% else 'equilibrium'). The model is open-loop unstable (growth rate a few 1/s):
% the pulse train makes it diverge, shorten tend to look at it.
vertical_uri = getenv('PLASMALESS_VERTICAL_EQUILIBRIUM_URI');
vertical_time = str2double(getenv('PLASMALESS_VERTICAL_EQUILIBRIUM_TIME'));   % NaN when unset
vertical_passive = getenv('PLASMALESS_VERTICAL_PASSIVE_CURRENTS');
if isempty(vertical_passive), vertical_passive = 'equilibrium'; end
if ~isempty(linear_model_dir) && ~isempty(vertical_uri)
    error('plasmaless:model','Set either the linear model directory or the vertical equilibrium URI, not both.');
end
if ~isempty(getenv('PLASMALESS_VERTICAL_EQUILIBRIUM_TIME')) && isnan(vertical_time)
    error('plasmaless:model','PLASMALESS_VERTICAL_EQUILIBRIUM_TIME is not a number.');
end
if isnan(vertical_time), vertical_time = []; end

%% load machine description
% DD4 entry shipped in data/md_dd4, made with tools/convert_md_dd3_to_dd4.py
% (IMAS-MATLAB DD-4.x module); also read for the linear model (coil names/order)
md_uri = ['imas:hdf5?path=' fullfile(repo_root,'data','md_dd4')];
fprintf('Loading machine description from %s\n', md_uri);
[em_coupling,pf_active_md,pf_passive_md] = load_plasmaless_md(md_uri);

%% build model and initial state (cold start)
if ~isempty(vertical_uri)
    % deviations from the equilibrium currents; the saved Ia, Ie and Y (Ip, rc,
    % zc) are absolute values, op + deviation (op.voltage = 0: the pulse train
    % is applied as is); do_plot uses plot_plasmaless_linear
    addpath(fullfile(repo_root,'vertical_instability'));
    eq_vertical = extract_equilibrium_data(vertical_uri, vertical_time, ...
        struct('passive_currents',vertical_passive,'pf_passive_md',pf_passive_md));
    [model, x] = configure_vertical_instability_model(dt, eq_vertical, em_coupling, ...
        pf_active_md, pf_passive_md);
    results_file = fullfile(repo_root,'plasmaless_vertical_results.mat');
elseif isempty(linear_model_dir)
    [model, x] = configure_plasmaless_model(em_coupling, pf_active_md, pf_passive_md, dt);
    results_file = fullfile(repo_root,'plasmaless_imas_results.mat');  % saved time, Ia, Ie, names
else
    % the pulse train is then a voltage deviation from the linearisation point
    % (op.voltage of the descriptor, 0 for the received model); the saved Ia
    % and Y are absolute values, op + deviation (op = 0 where not given)
    [model, x] = configure_linear_model(linear_model_dir, dt, pf_active_md);
    results_file = fullfile(repo_root,'plasmaless_linear_results.mat');
end
fprintf('Model built (%s): %d states, %d coils, %d passive loops, %d observables, dt = %g s\n', ...
    model.kind, model.n_states, model.n_coils, model.n_loops, numel(model.y_names), dt);
% warm start of the circuit model from a one-slice pf_active/pf_passive pair:
% [model, x] = configure_plasmaless_model(em_coupling, pf_active_md, pf_passive_md, dt, ...
%     pf_active_slice, pf_passive_slice);

%% prepare inputs
TS = timeseries_plasmaless_model(tstart,tend,dt);

%% run model (one IMAS step per time step, see simulator_plasmaless_model_imas.m)
if write_output, uri_arg = out_uri; else, uri_arg = []; end
[time, Ia, Ie, Y] = simulator_plasmaless_model_imas(TS, model, x, uri_arg);

%% save results and plot (plot_plasmaless_model.m / plot_plasmaless_linear.m)
coil_names = model.coil_names;
loop_names = model.loop_names;
y_names = model.y_names;      % empty for the circuit model
model_kind = model.kind;
save(results_file,'time','Ia','Ie','Y','coil_names','loop_names','y_names','model_kind','-v7');
fprintf('Results saved to %s\n', results_file);

% same output structure as the Simulink/MATLAB simulators of run_plasmaless_model.m:
% out.simout is a timeseries of [Ia Ie] (nt x (n_coils+n_loops)); TS, tstart, tend as above
out.tout = time;
out.simout = timeseries([Ia Ie], time);
if do_plot
    if strcmp(model_kind,'circuit')
        plot_plasmaless_model;
    else
        plot_plasmaless_linear;
    end
end
