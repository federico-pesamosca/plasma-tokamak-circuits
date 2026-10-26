function [pf_active_out, pf_passive_out, x_new, y_out] = plasmaless_step_ids(model, pf_active_in, x)
% one time step of the plasmaless model (as plasmaless_timestep.slx)
% model:        from configure_plasmaless_model (circuit), configure_linear_model
%               or configure_vertical_instability_model (vertical)
% pf_active_in: one slice (homogeneous_time 1) with the input voltages u at t:
%               the power-supply voltages (pf_active_in.supply, one per model
%               input, matched by short key against model.input_names) are
%               preferred; without supplies the coil voltages are used and
%               mapped with u = model.coil_to_input * v_coil (for ITER the VS1
%               input is then 0), see plasmaless_input_voltages
% x:            state at t (model.n_states x 1)
% Step: du = u - op.voltage; x_new = A*x + B*du; the outputs at t+dt are
%   coil currents  op.coil_current + x_to_coil*x_new   (pf_active_out.coil; the
%                  coil voltages are echoed only when the input carried coil
%                  voltages)
%   supply voltages u (pf_active_out.supply, names model.supply_names: the
%                  voltages applied over [t, t+dt])
%   loop currents  op.loop_current + x_to_loop*x_new   (pf_passive_out)
%   observables    op.y + C_y*x_new                    (y_out: names, values, known)
% For the circuit model op is zero and y_out is empty; for a linear model the
% states are deviations from the linearisation point op; for the vertical
% model they are deviations from the equilibrium currents (op.coil_current,
% op.loop_current) and y_out holds Ip, rc, zc.

if pf_active_in.ids_properties.homogeneous_time ~= 1 || numel(pf_active_in.time) ~= 1
    error('plasmaless:slice', ['pf_active input must have homogeneous_time 1 and ', ...
        'exactly one time slice (homogeneous_time %d, %d times).'], ...
        pf_active_in.ids_properties.homogeneous_time, numel(pf_active_in.time));
end
t = pf_active_in.time;
if ~isequal(size(x),[model.n_states 1])
    error('plasmaless:state','State must be a %dx1 vector.', model.n_states);
end

% input voltages in model (power-supply) order
[u, source] = plasmaless_input_voltages(model, pf_active_in);
if strcmp(source,'coil')
    % echo of the coil voltages of the input (same order as model.coil_names)
    names = cellfun(@plasmaless_coil_key, pf_active_in.coil(:)', 'UniformOutput', false);
    [~, idx] = ismember(model.coil_names, names);
    v_coil = cellfun(@(c) c.voltage.data, pf_active_in.coil(idx));
end

% state update: use the same function of simulator_plamaless_model (matlab/simulink)
% x_new = model.A*x + model.B*du;
du = u - model.op.voltage;
[~,x_new] = plasmaless_timestep_forward(du,x,model.A,model.B,model.C,model.D);
i_coil = model.op.coil_current + model.x_to_coil*x_new;
i_loop = model.op.loop_current + model.x_to_loop*x_new;
y_out.names = model.y_names;
y_out.values = model.op.y + model.C_y*x_new;
y_out.known = model.op.y_known;

% outputs
label = model_label(model);
pf_active_out = ids_init('pf_active');
pf_active_out.ids_properties.homogeneous_time = 1;
pf_active_out.ids_properties.comment = [label ': coil currents'];
pf_active_out.time = t + model.dt;
pf_active_out.coil = ids_allocate('pf_active','coil',model.n_coils);
for ii = 1:model.n_coils
    pf_active_out.coil{ii}.name = model.coil_names{ii};
    if isfield(pf_active_out.coil{ii},'identifier')
        pf_active_out.coil{ii}.identifier = model.coil_names{ii};
    end
    pf_active_out.coil{ii}.current.data = i_coil(ii);
    if strcmp(source,'coil')
        pf_active_out.coil{ii}.voltage.data = v_coil(ii);
    end
end
% the voltages applied over [t, t+dt], one per model input
pf_active_out.supply = ids_allocate('pf_active','supply',model.n_inputs);
for jj = 1:model.n_inputs
    pf_active_out.supply{jj}.name = model.supply_names{jj};
    if isfield(pf_active_out.supply{jj},'identifier')
        pf_active_out.supply{jj}.identifier = model.supply_names{jj};
    end
    pf_active_out.supply{jj}.voltage.data = u(jj);
end

pf_passive_out = ids_init('pf_passive');
pf_passive_out.ids_properties.homogeneous_time = 1;
pf_passive_out.ids_properties.comment = [label ': passive loop currents'];
pf_passive_out.time = t + model.dt;
pf_passive_out.loop = ids_allocate('pf_passive','loop',model.n_loops);
for ii = 1:model.n_loops
    pf_passive_out.loop{ii}.name = model.loop_names{ii};
    pf_passive_out.loop{ii}.current = i_loop(ii);
end

end

function label = model_label(model)
% text of the ids_properties.comment of the outputs
switch model.kind
    case 'circuit'
        label = 'plasmaless coil-vessel model';
    case 'linear'
        label = ['plasmaless linear model (' model.source.description ')'];
    case 'vertical'
        label = ['plasmaless vertical-instability model (' model.source.description ')'];
    otherwise
        error('plasmaless:model','Unknown model kind ''%s''.', model.kind);
end
end
