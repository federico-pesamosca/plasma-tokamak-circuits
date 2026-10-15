function [model, x0] = configure_plasmaless_model(em_coupling, pf_active, pf_passive, dt, ...
    pf_active_init, pf_passive_init)
% build the discrete plasmaless coil-vessel model from DD4 machine-description IDSs
% (see load_plasmaless_md.m, which reads the .h5 entry in data/md_dd4)
% em_coupling: DD4, coupling_matrix entries named
%              mutual_active_active, mutual_passive_active, mutual_passive_passive
% pf_active:   coils (name, resistance), same order as em_coupling
% pf_passive:  loops (name, resistance), same order as em_coupling
% The resistance of every coil and loop is read from the IDSs.
% Output model: kind 'circuit', sysd, A, B, C, D (discrete), dt, M, R, and the
% fields shared with configure_linear_model: n_states, state_names
% (= [coil_names loop_names]), coil_names, n_coils, x_to_coil (coil currents =
% x_to_coil*x), loop_names, n_loops, x_to_loop, y_names (none), C_y (0 x n),
% op (operating point, all zeros incl. loop_current: the circuit model is absolute).
% pf_active_init, pf_passive_init (optional, warm start): single time slice
% giving the initial state x0 = [Ia; Ie]; coils are matched to the model by
% short name, loops by position (name checked). Without them x0 is zero.

if nargin ~= 4 && nargin ~= 6
    error('plasmaless:nargin','Use 4 inputs (cold start) or 6 inputs (warm start).');
end

[Maa, Mea, Mee] = get_mutuals(em_coupling);
n_coils = numel(pf_active.coil);
n_loops = numel(pf_passive.loop);
if ~isequal(size(Maa),[n_coils n_coils]) || ~isequal(size(Mea),[n_loops n_coils]) ...
        || ~isequal(size(Mee),[n_loops n_loops])
    error('plasmaless:size', ['em_coupling sizes (Maa %s, Mea %s, Mee %s) ', ...
        'do not match %d pf_active coils and %d pf_passive loops.'], ...
        mat2str(size(Maa)), mat2str(size(Mea)), mat2str(size(Mee)), n_coils, n_loops);
end
M = [Maa Mea' ; Mea Mee];

% resistances
Raa = zeros(n_coils);
coil_names = cell(1,n_coils);
for ii = 1:n_coils
    coil_names{ii} = plasmaless_coil_key(pf_active.coil{ii});
    Raa(ii,ii) = get_resistance(pf_active.coil{ii}, ['pf_active coil ' coil_names{ii}]);
end
if numel(unique(coil_names)) ~= n_coils
    error('plasmaless:coilName','Duplicate pf_active coil names.');
end
Ree = zeros(n_loops);
loop_names = cell(1,n_loops);
for ii = 1:n_loops
    loop_names{ii} = strtrim(pf_passive.loop{ii}.name);
    Ree(ii,ii) = get_resistance(pf_passive.loop{ii}, ['pf_passive loop ' loop_names{ii}]);
end
R = blkdiag(Raa,Ree);

% build state space
Ac = -M\R; Bc = M\eye(n_coils+n_loops,n_coils);
Cc = eye(size(M)); Dc = zeros(size(Cc,1),size(Bc,2));
sys = ss(Ac,Bc,Cc,Dc);
sysd = c2d(sys,dt);

n_states = n_coils+n_loops;
model.kind = 'circuit';
model.sysd = sysd;
model.A = sysd.A;
model.B = sysd.B;
model.C = sysd.C;
model.D = sysd.D;
model.dt = dt;
model.M = M;
model.R = R;
model.n_states = n_states;
model.state_names = [coil_names loop_names];
model.coil_names = coil_names;
model.n_coils = n_coils;
model.x_to_coil = [eye(n_coils) zeros(n_coils,n_loops)];
model.loop_names = loop_names;
model.n_loops = n_loops;
model.x_to_loop = [zeros(n_loops,n_coils) eye(n_loops)];
model.y_names = {};
model.C_y = zeros(0,n_states);
% operating point (zero: states, inputs and outputs are absolute values)
model.op = struct('coil_current',zeros(n_coils,1),'coil_current_given',false(n_coils,1), ...
    'loop_current',zeros(n_loops,1),'voltage',zeros(n_coils,1),'voltage_given',false(n_coils,1),'y',zeros(0,1),'y_known',false(0,1));

% initial state
x0 = zeros(n_coils+n_loops,1);   % cold start, as x0 = zeros(size(A,1),1) in run_plasmaless_model.m
if nargin == 6
    idx = match_coils(model, pf_active_init);
    for ii = 1:n_coils
        x0(ii) = single_value(pf_active_init.coil{idx(ii)}.current.data, ...
            ['pf_active coil ' coil_names{ii} ' current']);
    end
    if numel(pf_passive_init.loop) ~= n_loops
        error('plasmaless:loops','pf_passive has %d loops, the model %d.', ...
            numel(pf_passive_init.loop), n_loops);
    end
    for ii = 1:n_loops
        name = strtrim(pf_passive_init.loop{ii}.name);
        if ~strcmp(name, loop_names{ii})
            error('plasmaless:loops','pf_passive loop %d is ''%s'', the model expects ''%s''.', ...
                ii, name, loop_names{ii});
        end
        x0(n_coils+ii) = single_value(pf_passive_init.loop{ii}.current, ...
            ['pf_passive loop ' name ' current']);
    end
end

end

function [Maa, Mea, Mee] = get_mutuals(em_coupling)
if ~isfield(em_coupling,'coupling_matrix')
    error('plasmaless:emCoupling','em_coupling has no coupling_matrix (DD4 expected).');
end
Maa = get_coupling_matrix(em_coupling,'mutual_active_active','pf_active/coil','pf_active/coil');
Mea = get_coupling_matrix(em_coupling,'mutual_passive_active','pf_passive/loop','pf_active/coil');
Mee = get_coupling_matrix(em_coupling,'mutual_passive_passive','pf_passive/loop','pf_passive/loop');
if isempty(Maa) || isempty(Mea) || isempty(Mee)
    error('plasmaless:emCoupling','em_coupling mutual inductance matrices are empty.');
end
end

function r = get_resistance(elem, label)
r = elem.resistance;
if ~isscalar(r) || r == -9e40 || ~isfinite(r)
    error('plasmaless:resistance','%s: resistance missing.', label);
end
end

function v = single_value(data, label)
if numel(data) ~= 1 || data == -9e40 || ~isfinite(data)
    error('plasmaless:data','%s: expected exactly one value, got %d.', label, numel(data));
end
v = data;
end

function idx = match_coils(model, pf_active)
% indices of the model coils in pf_active.coil, matched by short name
% error if a model coil is missing or pf_active has other coils

names = cellfun(@plasmaless_coil_key, pf_active.coil(:)', 'UniformOutput', false);
extra = setdiff(names, model.coil_names);
if ~isempty(extra) || numel(names) ~= model.n_coils
    error('plasmaless:coils','pf_active coils {%s} do not match the model coils {%s}.', ...
        strjoin(names,','), strjoin(model.coil_names,','));
end
[~, idx] = ismember(model.coil_names, names);
if any(idx == 0)
    error('plasmaless:coils','pf_active misses coils %s.', ...
        strjoin(model.coil_names(idx == 0),','));
end
end
