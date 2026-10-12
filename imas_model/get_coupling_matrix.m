function data = get_coupling_matrix(em_coupling, name, rows, cols)
% data of the DD4 em_coupling.coupling_matrix entry called name
% rows, cols: node that rows_uri/columns_uri must point at, e.g. 'pf_active/coil',
%             'pf_passive/loop' or 'grid_ggd(1)/space(1)/objects_per_dimension(1)/object';
%             accepted are 'node(:)' or the explicit list 'node(1)'..'node(n)' in order,
%             n = size(data,1) resp. size(data,2)
% The entry must be unique and of quantity magnetic_flux (index 1).
% Used by configure_plasmaless_model (conductor mutuals) and
% plasmaless_grid_from_ids (coil-to-grid flux maps).

if ~isfield(em_coupling,'coupling_matrix')
    error('plasmaless:emCoupling','em_coupling has no coupling_matrix (DD4 entry expected).');
end
found = [];
for ii = 1:numel(em_coupling.coupling_matrix)
    if strcmp(strtrim(em_coupling.coupling_matrix{ii}.name), name)
        found(end+1) = ii; %#ok<AGROW>
    end
end
if numel(found) ~= 1
    error('plasmaless:emCoupling', ['em_coupling has %d coupling_matrix named %s ', ...
        '(DD3 entries read with a DD4 API have none: convert them with ', ...
        'tools/convert_md_dd3_to_dd4.py).'], numel(found), name);
end
cm = em_coupling.coupling_matrix{found};
if cm.quantity.index ~= 1
    error('plasmaless:emCoupling','coupling_matrix %s: quantity index %d is not 1 (magnetic_flux).', ...
        name, cm.quantity.index);
end
data = cm.data;
if isempty(data)
    error('plasmaless:emCoupling','coupling_matrix %s: data is empty.', name);
end
check_uri(cm.rows_uri, rows, size(data,1), name);
check_uri(cm.columns_uri, cols, size(data,2), name);
end

function check_uri(uri, node, n, name)
% accept 'node(:)' or the explicit list 'node(1)'..'node(n)'
uri = cellstr(uri);
if isequal(uri, {[node '(:)']})
    return
end
expected = arrayfun(@(i) sprintf('%s(%d)',node,i), 1:n, 'UniformOutput', false);
if ~isequal(strtrim(uri(:)'), expected)
    error('plasmaless:emCoupling','coupling_matrix %s: rows/columns uri do not refer to %s(1..%d) in order.', ...
        name, node, n);
end
end
