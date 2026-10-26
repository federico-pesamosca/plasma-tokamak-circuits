function [model, x0] = configure_plasmaless_model(em_coupling, pf_active, pf_passive, dt, ...
    pf_active_init, pf_passive_init)
% build the discrete plasmaless coil-vessel model from DD4 machine-description IDSs
% (see load_plasmaless_md.m, which reads the .h5 entry in data/md_dd4)
% em_coupling: DD4, coupling_matrix entries named
%              mutual_active_active, mutual_passive_active, mutual_passive_passive
% pf_active:   coils (name, resistance), same order as em_coupling; its
%              supplies (name) give the supply names of the IMAS output
% pf_passive:  loops (name, resistance), same order as em_coupling
% The resistance of every coil and loop is read from the IDSs.
%
% The state is the current of the 12 active CIRCUITS and of the passive loops,
% x = [Ic; Ie], Ic = {CS3U CS2U CS1 CS2L CS3L PF1..PF6 VS3} (see
% plasmaless_circuit_topology): CS1U/CS1L are in series (Ia_CS1U = Ia_CS1L = Ic_CS1)
% and VS3U/VS3L in antiseries (Ia_VS3U = -Ic_VS3, Ia_VS3L = Ic_VS3). With Ia = T*Ic
% the circuit matrices are Mcc = T'*Maa*T, Mce = T'*Mae, Rcc = T'*Raa*T.
% The input is the 14 power-supply voltages Vps = {CS3U CS2U CS1U CS1L CS2L CS3L
% PF1..PF5 VS1 PF6 VS3} mapped to the circuit voltages Vc = P*Vps (see
% map_powersupplies_to_circuits); the output is y = [Ia (14 coils); Ie], i.e.
% the same 14 coil currents and passive loop currents as before.
% Output model: kind 'circuit', sysd, A, B, C, D (discrete), dt, M, R (circuits,
% [Mcc Mce; Mec Mee]), M_coils, R_coils (single coils), T, P, topology,
% coil_names, circuit_names, ps_names, loop_names, n_coils, n_circuits, n_ps,
% n_loops, n_states, and the fields shared with configure_linear_model and
% configure_vertical_instability_model (IMAS interface, plasmaless_step_ids):
% state_names (= [circuit_names loop_names]), x_to_coil (coil currents =
% x_to_coil*x = [T 0]*x), x_to_loop, n_inputs (= n_ps), input_names (= ps_names,
% column order of B), supply_names (pf_active supply names written in the
% output, see plasmaless_supply_names), coil_to_input (= topology.coil_to_ps,
% Vps from coil voltages when an input carries coil voltages only), y_names
% (none), C_y (0 x n_states), op (operating point, all zeros incl. loop_current;
% voltage/voltage_given sized n_inputs: the circuit model is absolute).
% pf_active_init, pf_passive_init (optional, warm start): single time slice
% giving the initial state x0 = [Ic; Ie]; coils are matched to the model by
% short name (the coil currents must satisfy the series/antiseries constraints),
% loops by position (name checked). Without them x0 is zero.

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

% circuits: coils in series/antiseries (CS1U+CS1L, VS3U-VS3L) share one current.
% Ia = T*Ic, state I = [Ic; Ie], full transformation Tf = blkdiag(T, I):
%   Mcc.. = Tf'*M*Tf,  R = Tf'*R*Tf,  Vc = T'*Va  (coil voltages) = P*Vps
topo = plasmaless_circuit_topology(coil_names);
n_circuits = topo.n_circuits;
Tf = blkdiag(topo.T, eye(n_loops));
Mc = Tf'*M*Tf;                 % [Mcc Mce; Mec Mee]
Rc = Tf'*R*Tf;                 % [Rcc 0; 0 Ree]
Mc = (Mc+Mc')/2;

% build state space: states [Ic; Ie], inputs Vps (power supplies),
% outputs [Ia (14 coils); Ie]
nx = n_circuits+n_loops;
S = eye(nx,n_circuits);        % circuit voltages enter the circuit rows only
Ac = -Mc\Rc;
Bc = Mc\(S*topo.P);
Cc = Tf;
Dc = zeros(size(Cc,1),size(Bc,2));
sys = ss(Ac,Bc,Cc,Dc);
sysd = c2d(sys,dt);

model.kind = 'circuit';
model.sysd = sysd;
model.A = sysd.A;
model.B = sysd.B;
model.C = sysd.C;
model.D = sysd.D;
model.dt = dt;
model.M = Mc;                  % circuit inductance [Mcc Mce; Mec Mee]
model.R = Rc;                  % circuit resistance [Rcc 0; 0 Ree]
model.M_coils = M;             % single-coil inductance [Maa Mae; Mea Mee]
model.R_coils = R;
model.coil_names = coil_names;
model.loop_names = loop_names;
model.circuit_names = topo.circuit_names;
model.ps_names = topo.ps_names;
model.topology = topo;
model.T = topo.T;              % Ia = T*Ic
model.P = topo.P;              % Vc = P*Vps
model.n_coils = n_coils;
model.n_circuits = n_circuits;
model.n_ps = topo.n_ps;
model.n_loops = n_loops;
model.n_states = nx;
% IMAS interface (plasmaless_step_ids): states, inputs and outputs by name
model.state_names = [topo.circuit_names loop_names];
model.x_to_coil = [topo.T zeros(n_coils,n_loops)];         % Ia = x_to_coil*x
model.x_to_loop = [zeros(n_loops,n_circuits) eye(n_loops)];
model.n_inputs = topo.n_ps;
model.input_names = topo.ps_names;                         % column order of B
model.supply_names = plasmaless_supply_names(model.input_names, pf_active);
model.coil_to_input = topo.coil_to_ps;                     % Vps = coil_to_input*Va
model.y_names = {};
model.C_y = zeros(0,nx);
% operating point (zero: states, inputs and outputs are absolute values)
model.op = struct('coil_current',zeros(n_coils,1),'coil_current_given',false(n_coils,1), ...
    'loop_current',zeros(n_loops,1),'voltage',zeros(model.n_inputs,1), ...
    'voltage_given',false(model.n_inputs,1),'y',zeros(0,1),'y_known',false(0,1));

% initial state [Ic; Ie] (circuit currents)
x0 = zeros(nx,1);   % cold start
if nargin == 6
    idx = match_coils(model, pf_active_init);
    Ia0 = zeros(1,n_coils);
    for ii = 1:n_coils
        Ia0(ii) = single_value(pf_active_init.coil{idx(ii)}.current.data, ...
            ['pf_active coil ' coil_names{ii} ' current']);
    end
    x0(1:n_circuits) = map_coils_to_circuits(Ia0,topo);
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
        x0(n_circuits+ii) = single_value(pf_passive_init.loop{ii}.current, ...
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
