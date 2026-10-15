% run simulator_plasmaless.slx (simulink) 
% or simulator_plasmaless_model_matlab.m (matlab)
% featuring dynamical plasmaless coil-vessel circuit model
% (or an external linear model, see linear_model_dir, or the vertical-instability
% model, see vertical_uri) and plot the results

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
% '' = circuit model built from the machine description; the path of a directory
% holding linear_model.json + matrices = external linear model (deviations from
% its linearisation point), e.g. fullfile(fileparts(mfilename('fullpath')),
% 'linearModelSimple_105073_t84.417598'). Default: the environment variable
% PLASMALESS_LINEAR_MODEL_DIR when set, else ''.
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

%% load model (DD4 machine description in data/md_dd4)
[em_coupling,pf_active,pf_passive] = load_plasmaless_md();
if ~isempty(vertical_uri)
    % the outputs are deviations from the equilibrium currents (no op offsets)
    eq_vertical = extract_equilibrium_data(vertical_uri, vertical_time, ...
        struct('passive_currents',vertical_passive,'pf_passive_md',pf_passive));
    [model,x0] = configure_vertical_instability_model(dt,eq_vertical,em_coupling,pf_active,pf_passive);
    results_file = fullfile(fileparts(mfilename('fullpath')),'plasmaless_matlab_vertical_results.mat');
elseif isempty(linear_model_dir)
    [model,x0] = configure_plasmaless_model(em_coupling,pf_active,pf_passive,dt);
    results_file = fullfile(fileparts(mfilename('fullpath')),'plasmaless_matlab_results.mat');
else
    [model,x0] = configure_linear_model(linear_model_dir,dt,pf_active);
    results_file = fullfile(fileparts(mfilename('fullpath')),'plasmaless_matlab_linear_results.mat');
end
% matrices used by the simulink model and by the matlab simulator (base workspace):
% the output is [coil currents; loop currents; observables] (for the circuit model
% C_full = eye(n_states), as before). Both simulators feed u directly, without the
% op.voltage subtraction of plasmaless_step_ids: op is zero for the circuit model;
% for a linear model the pulse train is a voltage deviation and the outputs are
% deviations too (no op offsets added).
C_full = [model.x_to_coil; model.x_to_loop; model.C_y];
D_full = zeros(size(C_full,1), model.n_coils);

%% prepare timeseries
TS = timeseries_plasmaless_model(tstart,tend,dt);

%% run model
switch mode
    case 'simulink'
        A = model.A; B = model.B; C = C_full; D = D_full;
        out = sim('simulator_plasmaless_model_simulink.slx');
    case 'matlab'
        model_sim = model; model_sim.C = C_full; model_sim.D = D_full;
        out = simulator_plasmaless_model_matlab(TS,model_sim,x0);
end

%% save results (same variables and units as run_plasmaless_imas.m)
time = out.simout.Time;
y = squeeze(out.simout.Data);
if size(y,1) ~= numel(time), y = y.'; end
n_coils = model.n_coils;
n_loops = model.n_loops;
Ia = y(:,1:n_coils);
Ie = y(:,n_coils+(1:n_loops));
Y = y(:,n_coils+n_loops+(1:numel(model.y_names)));
coil_names = model.coil_names;
loop_names = model.loop_names;
y_names = model.y_names;      % empty for the circuit model
model_kind = model.kind;
save(results_file,'time','Ia','Ie','Y','coil_names','loop_names','y_names','model_kind','-v7');
fprintf('Results saved to %s\n', results_file);

%% plot model
if do_plot
    if strcmp(model_kind,'circuit')
        plot_plasmaless_model;
    else
        plot_plasmaless_linear;
    end
end
