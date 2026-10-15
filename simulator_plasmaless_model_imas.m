function [time, Ia, Ie, Y] = simulator_plasmaless_model_imas(TS, model, x0, out_uri)
% run the plasmaless model step by step through its IMAS interface
% TS:    timeseries of the coil voltages (nt x n_coils), see timeseries_plasmaless_model.m
% model: model from configure_plasmaless_model, configure_linear_model or
%        configure_vertical_instability_model
% x0:    initial state (model.n_states x 1)
% out_uri (optional): IMAS URI where the pf_active/pf_passive output slices are
%        written (ids_put, then ids_put_slice); nothing is written if omitted or empty
% time:  TS.Time (nt x 1)
% Ia, Ie: coil (nt x n_coils) and passive loop (nt x n_loops) currents, row k at time(k)
% Y:     observables (nt x numel(model.y_names), nt x 0 for the circuit model);
%        row 1 is op + the output matrices applied to x0, later rows the step outputs

write_output = nargin > 3 && ~isempty(out_uri);

%% prepare inputs
pf_active_in = ids_init('pf_active');
pf_active_in.ids_properties.homogeneous_time = 1;
pf_active_in.coil = ids_allocate('pf_active','coil',model.n_coils);
for ii = 1:model.n_coils
    pf_active_in.coil{ii}.name = model.coil_names{ii};
end

%% run model
x = x0;
time = TS.Time;
nt = numel(time);
Ia = zeros(nt,model.n_coils);
Ie = zeros(nt,model.n_loops);
Y = zeros(nt,numel(model.y_names));
Ia(1,:) = (model.op.coil_current + model.x_to_coil*x)';
Ie(1,:) = (model.op.loop_current + model.x_to_loop*x)';
Y(1,:) = (model.op.y + model.C_y*x)';
if write_output
    idx_out = imas_open(out_uri,43);
end
fprintf('Running %d steps, t = %g .. %g s\n', nt-1, time(1), time(end));
tic;
for k = 1:nt-1
    % one-slice input, as sent by a coupled controller
    pf_active_in.time = time(k);
    for ii = 1:model.n_coils
        pf_active_in.coil{ii}.voltage.data = TS.Data(k,ii);
    end
    [pf_active_out, pf_passive_out, x, y_out] = plasmaless_step_ids(model, pf_active_in, x);
    Ia(k+1,:) = cellfun(@(c) c.current.data, pf_active_out.coil);
    Ie(k+1,:) = cellfun(@(l) l.current, pf_passive_out.loop);
    Y(k+1,:) = y_out.values';
    if floor(time(k+1)+1e-9) > floor(time(k)+1e-9)  % progress every simulated second
        fprintf('t = %6.3f s (%3.0f%%), elapsed %.1f s\n', time(k+1), 100*k/(nt-1), toc);
        drawnow;
    end
    if write_output
        if k == 1
            ids_put(idx_out,'pf_active',pf_active_out);
            ids_put(idx_out,'pf_passive',pf_passive_out);
        else
            ids_put_slice(idx_out,'pf_active',pf_active_out);
            ids_put_slice(idx_out,'pf_passive',pf_passive_out);
        end
    end
end
if write_output
    imas_close(idx_out);
end
fprintf('Simulation done in %.1f s wall time\n', toc);

end
