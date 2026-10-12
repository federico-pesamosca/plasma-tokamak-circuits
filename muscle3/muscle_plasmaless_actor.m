function muscle_plasmaless_actor()
% MUSCLE3 actor around the plasmaless coil+vessel model (imas_model/), a
% drop-in replacement of NICE direct evolutive (nice_evo_rd) for controller
% tests in PDS (workflow plasmaless_controller). Run with
%   matlab -batch "pyenv('Version',<python with muscle3 0.10>,'ExecutionMode','InProcess'); ...
%                  addpath(<repo>/muscle3); muscle_plasmaless_actor"
% See muscle3/README.md for ports, settings, timing and assumptions.
%
% Ports
%   F_INIT equilibrium_f_init  multi-time equilibrium (reference trajectory);
%                              its message timestamp is the start time t0
%          pf_active_f_init    multi-time pf_active (machine description +
%                              scenario coil currents/voltages)
%   S      pf_active_s         one-slice pf_active with the coil voltages
%                              (from the controller, coil names absent)
%   O_I    equilibrium_o_i     one-slice reference pass-through equilibrium
%          pf_active_o_i       one-slice coil currents (model)
%          pf_passive_o_i      one-slice passive loop currents (model)
% Settings: dt, t_interval (must equal dt), t_end (mandatory); md_uri
% (optional, default: data/md_dd4 of this repository).
%
% Timing mirrors nice_imas_evo_rd_muscle3 (NICE 3.0.0.dev446,
% main_imas_evo_rd_muscle3.cc): t0 = timestamp of equilibrium_f_init (cc:146),
% F_INIT pf_active sliced at t0 with the closest sample (cc:156); the first
% step uses that slice; after every step the outputs at t are sent with
% next_timestamp t+t_interval, or None when t+t_interval > t_end+1e-9, after
% which the actor stops (cc:178-195); otherwise pf_active_s is received and
% gives the voltages of the next step (cc:197).

% repository root (plasmaless_timestep_forward.m, the step function shared by
% all model variants, lives there) and imas_model/ next to this folder (not a
% setting)
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repo_root);
addpath(fullfile(repo_root,'imas_model'));

ports = py.dict();
ports{py.getattr(@py.ymmsl.Operator,"F_INIT")} = py.list({"equilibrium_f_init", "pf_active_f_init"});
ports{py.getattr(@py.ymmsl.Operator,"S")} = py.list({"pf_active_s"});
ports{py.getattr(@py.ymmsl.Operator,"O_I")} = py.list({"equilibrium_o_i", "pf_active_o_i", "pf_passive_o_i"});
flags = py.getattr(@py.libmuscle.InstanceFlags,"KEEPS_NO_STATE_FOR_NEXT_USE");
instance = py.libmuscle.Instance(ports, flags);

try
    while instance.reuse_instance()
        run_one(instance);
    end
catch ME
    fprintf(2,'muscle_plasmaless_actor: %s\n', getReport(ME,'extended','hyperlinks','off'));
    instance.error_shutdown(ME.message);
    rethrow(ME);
end
fprintf('muscle_plasmaless_actor: done\n');

end

function run_one(instance)
% one reuse of the instance: F_INIT, model build, time loop

% settings
dt = double(instance.get_setting('dt'));
t_interval = double(instance.get_setting('t_interval'));
t_end = double(instance.get_setting('t_end'));
% default machine description: data/md_dd4 shipped in this repository
md_uri = char(instance.get_setting('md_uri', pyargs('default', ['imas:hdf5?path=' ...
    fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','md_dd4')])));
if ~(dt > 0)
    error('plasmaless_actor:settings','dt must be > 0 (got %g).', dt);
end
if t_interval ~= dt
    % the model output is only defined at t0 + k*dt; an exchange period other
    % than dt would need sub-stepping that NICE does but this actor does not
    error('plasmaless_actor:settings','t_interval (%g) must equal dt (%g).', t_interval, dt);
end

% F_INIT
msg_eq = instance.receive("equilibrium_f_init");
t0 = double(msg_eq.timestamp);
eq_ref = imas_deserialize(uint8(msg_eq.data),'equilibrium');
msg_pfa = instance.receive("pf_active_f_init");
pfa_init = imas_deserialize(uint8(msg_pfa.data),'pf_active');
fprintf('muscle_plasmaless_actor: F_INIT received, t0 = %.10g, t_end = %.10g, dt = %g\n', t0, t_end, dt);

[eq_time, ip_ref, r_ref, z_ref] = reference_trajectory(eq_ref);
if ~(t_end > t0)
    error('plasmaless_actor:t_end','t_end (%.10g) must be > t0 (%.10g).', t_end, t0);
end
if t0 < eq_time(1) || t_end > eq_time(end)
    error('plasmaless_actor:t_end', ['[t0, t_end] = [%.10g, %.10g] is not inside the ', ...
        'F_INIT equilibrium time range [%.10g, %.10g].'], t0, t_end, eq_time(1), eq_time(end));
end

% machine description and model
fprintf('muscle_plasmaless_actor: loading machine description from %s\n', md_uri);
idx = imas_open(md_uri,40);
em_coupling = ids_get(idx,'em_coupling');
pf_active_md = ids_get(idx,'pf_active');
pf_passive_md = ids_get(idx,'pf_passive');
imas_close(idx);

% The controller reads and writes the coils by index (KCURR_RZIp
% muscle_IDS_NICE_output.m:94-97, muscle_NICE_input.m:45-51), in the order of
% the F_INIT pf_active; the model outputs are in machine-description order.
% Both must be the same coils in the same order (same short names), otherwise
% currents go to the wrong coils; any difference is an error.
names_init = cellfun(@short_name, pfa_init.coil(:)', 'UniformOutput', false);
if ~isequal(names_init, model_coil_names(pf_active_md))
    error('plasmaless_actor:coilOrder', ['F_INIT pf_active coils {%s} differ from the ', ...
        'machine-description coils {%s} (names or order).'], ...
        strjoin(names_init,','), strjoin(model_coil_names(pf_active_md),','));
end
coil_names = model_coil_names(pf_active_md);

pfa_t0 = closest_slice(pfa_init, t0, coil_names);
% ASSUMPTION: the vessel (pf_passive loop) currents are zero at t0. The
% scenario data carried by the workflow has no passive-structure currents, so
% there is no source for them; this is a modelling assumption, not a default
% of the model (configure_plasmaless_model needs them explicitly).
pfp_t0 = ids_init('pf_passive');
pfp_t0.ids_properties.homogeneous_time = 1;
pfp_t0.time = t0;
pfp_t0.loop = ids_allocate('pf_passive','loop',numel(pf_passive_md.loop));
for ii = 1:numel(pf_passive_md.loop)
    pfp_t0.loop{ii}.name = pf_passive_md.loop{ii}.name;
    pfp_t0.loop{ii}.current = 0;
end
[model, x] = configure_plasmaless_model(em_coupling, pf_active_md, pf_passive_md, dt, pfa_t0, pfp_t0);
fprintf('muscle_plasmaless_actor: model built, %d coils, %d passive loops\n', model.n_coils, model.n_loops);

% time loop
pfa_in = pfa_t0;           % first step: scenario voltages at t0 (as NICE)
t = t0;
n_sent = 0;
while true
    [pfa_out, pfp_out, x] = plasmaless_step_ids(model, pfa_in, x);
    t = pfa_out.time;      % = t + dt
    eq_out = reference_equilibrium(t, eq_time, ip_ref, r_ref, z_ref);

    more_after_this = t + t_interval <= t_end + 1e-9;   % as cc:184
    if more_after_this
        t_next = py.float(t + t_interval);
    else
        t_next = py.None;
    end
    send_ids(instance, "equilibrium_o_i", eq_out, 'equilibrium', t, t_next);
    send_ids(instance, "pf_active_o_i", pfa_out, 'pf_active', t, t_next);
    send_ids(instance, "pf_passive_o_i", pfp_out, 'pf_passive', t, t_next);
    n_sent = n_sent + 1;
    if ~more_after_this
        break
    end

    % next input; the stop decision is ours (t_end), not taken from this
    % message (the controller always sends next_timestamp None,
    % muscle_NICE_input.m:57)
    msg = instance.receive("pf_active_s");
    pfa_in = imas_deserialize(uint8(msg.data),'pf_active');
    t_msg = double(msg.timestamp);
    % The controller stamps its command with the time of the measurement it
    % answers (muscle_NICE_input.m:44,57); a different time means the loop is
    % out of step. 1e-6 s tolerance: well below dt, above round-off of the
    % accumulated time in Simulink/here.
    if abs(t_msg - t) > 1e-6 || numel(pfa_in.time) ~= 1 || abs(pfa_in.time - t) > 1e-6
        error('plasmaless_actor:time', ['pf_active_s at timestamp %.10g (IDS time %s) does not ', ...
            'answer the measurement sent at t = %.10g.'], t_msg, mat2str(pfa_in.time), t);
    end
    pfa_in.time = t;       % keep this actor's own clock (differences < 1e-6 s)
    % The controller does not send coil names; set them by index (F_INIT
    % order, checked above to be the machine-description order), as
    % nice_imas_evo_rd_muscle3 copies them by index from its own pf_active
    % (cc:206-209, "the controller doesn't give the name of the coils")
    if numel(pfa_in.coil) ~= numel(coil_names)
        error('plasmaless_actor:coils','pf_active_s has %d coils, the model %d.', ...
            numel(pfa_in.coil), numel(coil_names));
    end
    for ii = 1:numel(pfa_in.coil)
        pfa_in.coil{ii}.name = coil_names{ii};
    end
end
fprintf('muscle_plasmaless_actor: %d exchanges, last t = %.10g\n', n_sent, t);

end

function [eq_time, ip_ref, r_ref, z_ref] = reference_trajectory(eq)
% ip and boundary geometric axis of every F_INIT equilibrium time slice
if eq.ids_properties.homogeneous_time ~= 1
    error('plasmaless_actor:equilibrium','F_INIT equilibrium must have homogeneous_time 1 (got %d).', ...
        eq.ids_properties.homogeneous_time);
end
eq_time = eq.time(:);
n = numel(eq.time_slice);
if n ~= numel(eq_time) || n == 0
    error('plasmaless_actor:equilibrium','F_INIT equilibrium has %d time slices and %d times.', ...
        n, numel(eq_time));
end
if any(diff(eq_time) <= 0)
    error('plasmaless_actor:equilibrium','F_INIT equilibrium times are not increasing.');
end
ip_ref = zeros(n,1); r_ref = zeros(n,1); z_ref = zeros(n,1);
for ii = 1:n
    s = eq.time_slice{ii};
    ip_ref(ii) = check_value(s.global_quantities.ip, 'global_quantities.ip', eq_time(ii));
    r_ref(ii) = check_value(s.boundary.geometric_axis.r, 'boundary.geometric_axis.r', eq_time(ii));
    z_ref(ii) = check_value(s.boundary.geometric_axis.z, 'boundary.geometric_axis.z', eq_time(ii));
end
end

function v = check_value(v, label, t)
if numel(v) ~= 1 || v == -9e40 || ~isfinite(v)
    error('plasmaless_actor:equilibrium','F_INIT equilibrium %s missing at t = %.10g.', label, t);
end
end

function eq = reference_equilibrium(t, eq_time, ip_ref, r_ref, z_ref)
% one-slice equilibrium: ip and geometric axis of the reference trajectory,
% linearly interpolated at t. There is no plasma in this model; these values
% make the controller's Ip/R/Z errors zero (reference pass-through).
if numel(eq_time) == 1
    if t ~= eq_time
        error('plasmaless_actor:equilibrium','t = %.10g outside the one-slice F_INIT equilibrium.', t);
    end
    vals = [ip_ref r_ref z_ref];
else
    vals = interp1(eq_time, [ip_ref r_ref z_ref], t, 'linear');
end
if any(~isfinite(vals))
    error('plasmaless_actor:equilibrium','t = %.10g outside the F_INIT equilibrium time range.', t);
end
eq = ids_init('equilibrium');
eq.ids_properties.homogeneous_time = 1;
eq.ids_properties.comment = 'plasmaless model: no plasma, reference pass-through of ip and geometric axis';
eq.time = t;
eq.time_slice = ids_allocate('equilibrium','time_slice',1);
eq.time_slice{1}.time = t;
eq.time_slice{1}.global_quantities.ip = vals(1);
eq.time_slice{1}.boundary.geometric_axis.r = vals(2);
eq.time_slice{1}.boundary.geometric_axis.z = vals(3);
eq.code.output_flag = 0;
end

function names = model_coil_names(pf_active_md)
% short coil names of the machine description, in its order (the model order)
names = cellfun(@short_name, pf_active_md.coil(:)', 'UniformOutput', false);
end

function pfa = closest_slice(pfa_all, t, coil_names)
% one-slice pf_active at the sample of pfa_all closest to t (as NICE's
% getSlice(..., CLOSEST_SAMPLE), cc:156; equidistant samples: the earlier one)
% with coil currents and voltages, coils named coil_names (by index)
if pfa_all.ids_properties.homogeneous_time ~= 1
    error('plasmaless_actor:pf_active','F_INIT pf_active must have homogeneous_time 1 (got %d).', ...
        pfa_all.ids_properties.homogeneous_time);
end
tt = pfa_all.time(:);
if isempty(tt)
    error('plasmaless_actor:pf_active','F_INIT pf_active has no time.');
end
[~, k] = min(abs(tt - t));
pfa = ids_init('pf_active');
pfa.ids_properties.homogeneous_time = 1;
pfa.time = t;
n = numel(pfa_all.coil);
pfa.coil = ids_allocate('pf_active','coil',n);
for ii = 1:n
    c = pfa_all.coil{ii};
    if numel(c.current.data) ~= numel(tt) || numel(c.voltage.data) ~= numel(tt)
        error('plasmaless_actor:pf_active', ['F_INIT pf_active coil %d: %d currents and %d ', ...
            'voltages for %d times.'], ii, numel(c.current.data), numel(c.voltage.data), numel(tt));
    end
    pfa.coil{ii}.name = coil_names{ii};
    pfa.coil{ii}.current.data = c.current.data(k);
    pfa.coil{ii}.voltage.data = c.voltage.data(k);
end
fprintf('muscle_plasmaless_actor: initial coil currents/voltages from F_INIT pf_active sample t = %.10g\n', tt(k));
end

function send_ids(instance, port, ids, ids_name, t, t_next)
% serialization as KCURR_RZIp/muscle_NICE_input.m:53-57
serialized = imas_serialize(ids, ids_name);
s = sprintf('%d ' ,uint8(serialized));
x = py.numpy.fromstring(s, py.numpy.int8, int8(-1), char(' '));
serialized = py.bytes(x);
msg = py.libmuscle.Message(t, t_next, serialized);
instance.send(port, msg);
end

function key = short_name(coil)
% short coil key, same rule as imas_model/configure_plasmaless_model.m coil_key
if isfield(coil,'identifier') && ~isempty(strtrim(coil.identifier))
    key = coil.identifier;
else
    key = coil.name;
    tok = regexp(key,'\(([^()]*)\)','tokens','once');
    if ~isempty(tok)
        key = tok{1};
    end
end
key = strtrim(key);
end
