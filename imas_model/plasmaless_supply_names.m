function names = plasmaless_supply_names(input_names, pf_active_md)
% pf_active.supply names to write in the output IDS for the model inputs
% input_names:  1 x n_inputs cellstr of input (power-supply) short keys
% pf_active_md: machine-description pf_active; its supply{j}.name are matched
%               to input_names by short key (plasmaless_supply_key)
% names: 1 x n_inputs cellstr, pf_active_md.supply{j}.name of the supply whose
%        key is input_names{k}; where the machine description has no supplies
%        or no supply with that key, the input key itself.
% Shared by configure_plasmaless_model, configure_linear_model (and, through
% the base model, configure_vertical_instability_model).

names = input_names(:)';
if ~isfield(pf_active_md,'supply') || ~iscell(pf_active_md.supply) || isempty(pf_active_md.supply)
    return
end
md_names = cellfun(@(s) strtrim(s.name), pf_active_md.supply(:)', 'UniformOutput', false);
md_keys = cellfun(@plasmaless_supply_key, pf_active_md.supply(:)', 'UniformOutput', false);
for k = 1:numel(names)
    j = find(strcmp(md_keys, names{k}), 1);
    if ~isempty(j)
        names{k} = md_names{j};
    end
end

end
