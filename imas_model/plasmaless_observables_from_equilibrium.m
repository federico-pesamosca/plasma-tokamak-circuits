function [values, known] = plasmaless_observables_from_equilibrium(slice, names)
% read the observables names (vocabulary of plasmaless_observables_to_equilibrium)
% from one equilibrium time_slice struct (DD4 fields, see that function)
% values: numel(names) x 1, NaN where unknown; known: logical, true when the
% field exists, is non-empty, finite and not -9e40 (gap<k>: numel(gap) >= k;
% rX/zX: the first contour_tree node with critical_type 1; Dist2ndSep: never,
% DD4 has no field for it)

names = names(:)';
values = nan(numel(names),1);
known = false(numel(names),1);
for k = 1:numel(names)
    switch names{k}
        case 'rmag'
            v = get_path(slice,{'global_quantities','magnetic_axis','r'});
        case 'zmag'
            v = get_path(slice,{'global_quantities','magnetic_axis','z'});
        case 'rX'
            v = x_point(slice,'r');
        case 'zX'
            v = x_point(slice,'z');
        case 'Ip'
            v = get_path(slice,{'global_quantities','ip'});
        case 'Dist2ndSep'
            v = [];
        otherwise
            tok = regexp(names{k},'^gap([1-9][0-9]*)$','tokens','once');
            if isempty(tok)
                error('plasmaless:observable','Unknown observable ''%s''.', names{k});
            end
            kg = str2double(tok{1});
            v = [];
            gap = get_path(slice,{'boundary','gap'});
            if iscell(gap) && numel(gap) >= kg
                v = get_path(gap{kg},{'value'});
            end
    end
    if valid_scalar(v)
        values(k) = v;
        known(k) = true;
    end
end

end

function v = get_path(s, path)
% nested field value, [] when a field is missing
v = s;
for ii = 1:numel(path)
    if ~isstruct(v) || ~isfield(v,path{ii})
        v = [];
        return
    end
    v = v.(path{ii});
end
end

function v = x_point(slice, field)
% r or z of the first contour_tree node with critical_type 1 (saddle)
v = [];
nodes = get_path(slice,{'contour_tree','node'});
if ~iscell(nodes)
    return
end
for ii = 1:numel(nodes)
    ct = get_path(nodes{ii},{'critical_type'});
    if valid_scalar(ct) && ct == 1
        v = get_path(nodes{ii},{field});
        return
    end
end
end

function ok = valid_scalar(v)
ok = isnumeric(v) && isscalar(v) && isfinite(v) && v ~= -9e40;
end
