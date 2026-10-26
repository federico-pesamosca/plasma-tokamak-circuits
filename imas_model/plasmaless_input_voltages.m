function [u, source, given] = plasmaless_input_voltages(model, pf_active_in)
% input voltages u (model.n_inputs x 1, power-supply space) of one time step
% from a one-slice pf_active (no time checks here, see plasmaless_step_ids)
% model:        from configure_plasmaless_model, configure_linear_model or
%               configure_vertical_instability_model (input_names, coil_names,
%               coil_to_input)
% pf_active_in: one slice; supply voltages are preferred: when pf_active_in.supply
%               is non-empty every model input must be present exactly once
%               (matched by plasmaless_supply_key against model.input_names) and
%               u(k) = supply voltage k; otherwise the coil voltages are used
%               (matched by plasmaless_coil_key against model.coil_names) and
%               u = model.coil_to_input * v_coil (for ITER the VS1 input is then
%               0: the star-point converter voltage is not expressible from the
%               coil voltages).
% source: 'supply' or 'coil'
% given:  n_inputs x 1 logical, true where u(k) is determined by the slice
%         ('supply': all; 'coil': the inputs with a non-zero coil_to_input row)
% Errors: 'plasmaless:supplies' (missing, extra or duplicate supply),
%         'plasmaless:coils' (missing or extra coil), 'plasmaless:voltage'
%         (a voltage that is not exactly one finite value).

t = pf_active_in.time;
has_supplies = isfield(pf_active_in,'supply') && iscell(pf_active_in.supply) ...
    && numel(pf_active_in.supply) > 0;
if has_supplies
    idx = match_supplies(model, pf_active_in);
    u = zeros(model.n_inputs,1);
    for k = 1:model.n_inputs
        v = pf_active_in.supply{idx(k)}.voltage.data;
        if numel(v) ~= 1 || v == -9e40 || ~isfinite(v)
            error('plasmaless:voltage','pf_active supply %s: expected one voltage value at t = %g, got %d.', ...
                model.input_names{k}, t, numel(v));
        end
        u(k) = v;
    end
    source = 'supply';
    given = true(model.n_inputs,1);
else
    idx = match_coils(model, pf_active_in);
    v_coil = zeros(model.n_coils,1);
    for ii = 1:model.n_coils
        v = pf_active_in.coil{idx(ii)}.voltage.data;
        if numel(v) ~= 1 || v == -9e40 || ~isfinite(v)
            error('plasmaless:voltage','pf_active coil %s: expected one voltage value at t = %g, got %d.', ...
                model.coil_names{ii}, t, numel(v));
        end
        v_coil(ii) = v;
    end
    u = model.coil_to_input*v_coil;
    source = 'coil';
    given = any(model.coil_to_input ~= 0, 2);
end

end

function idx = match_supplies(model, pf_active)
% indices of the model inputs in pf_active.supply, matched by short key
% error if a model input is missing or pf_active has other/duplicate supplies

keys = cellfun(@plasmaless_supply_key, pf_active.supply(:)', 'UniformOutput', false);
extra = setdiff(keys, model.input_names);
if ~isempty(extra) || numel(keys) ~= model.n_inputs || numel(unique(keys)) ~= numel(keys)
    error('plasmaless:supplies','pf_active supplies {%s} do not match the model inputs {%s}.', ...
        strjoin(keys,','), strjoin(model.input_names,','));
end
[~, idx] = ismember(model.input_names, keys);
if any(idx == 0)
    error('plasmaless:supplies','pf_active misses supplies %s (inputs {%s}).', ...
        strjoin(model.input_names(idx == 0),','), strjoin(model.input_names,','));
end
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
