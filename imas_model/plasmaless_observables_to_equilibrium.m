function eq = plasmaless_observables_to_equilibrium(eq, y)
% write the known observables of a linear model into a one-slice equilibrium IDS
% y: struct names/values/known as returned by plasmaless_step_ids (absolute
%    values); entries with known false are not written
% DD4 (4.1.x) fields of time_slice{1}:
%   rmag, zmag  global_quantities.magnetic_axis.r / .z
%   rX, zX      contour_tree.node{1}.r / .z with critical_type 1 (saddle =
%               X-point): DD4 has no boundary_separatrix/x_point any more
%   Ip          global_quantities.ip
%   gap<k>      boundary.gap{k}.value, name 'gap<k>' (gap allocated to the largest k)
%   Dist2ndSep  not written: DD4 removed boundary_secondary_separatrix/
%               distance_inner_outer without a replacement
% The inverse read is plasmaless_observables_from_equilibrium.

if numel(eq.time_slice) ~= 1
    error('plasmaless:slice','equilibrium must have exactly one time slice (got %d).', numel(eq.time_slice));
end
names = y.names(:)';
known = logical(y.known(:)');
values = y.values(:)';
if numel(known) ~= numel(names) || numel(values) ~= numel(names)
    error('plasmaless:observable','y.names, y.values and y.known must have the same length.');
end
s = eq.time_slice{1};

% gaps: allocate once to the largest known k
kmax = 0;
for k = find(known)
    tok = regexp(names{k},'^gap([1-9][0-9]*)$','tokens','once');
    if ~isempty(tok)
        kmax = max(kmax,str2double(tok{1}));
    end
end
if kmax > numel(s.boundary.gap)
    gap = ids_allocate('equilibrium','time_slice/boundary/gap',kmax);
    gap(1:numel(s.boundary.gap)) = s.boundary.gap;
    s.boundary.gap = gap;
end

for k = find(known)
    v = values(k);
    switch names{k}
        case 'rmag'
            s.global_quantities.magnetic_axis.r = v;
        case 'zmag'
            s.global_quantities.magnetic_axis.z = v;
        case {'rX','zX'}
            if isempty(s.contour_tree.node)
                s.contour_tree.node = ids_allocate('equilibrium','time_slice/contour_tree/node',1);
            end
            s.contour_tree.node{1}.critical_type = int32(1);
            if strcmp(names{k},'rX')
                s.contour_tree.node{1}.r = v;
            else
                s.contour_tree.node{1}.z = v;
            end
        case 'Ip'
            s.global_quantities.ip = v;
        case 'Dist2ndSep'
            % no DD4 field (see header)
        otherwise
            tok = regexp(names{k},'^gap([1-9][0-9]*)$','tokens','once');
            if isempty(tok)
                error('plasmaless:observable','Unknown observable ''%s''.', names{k});
            end
            kg = str2double(tok{1});
            s.boundary.gap{kg}.name = names{k};
            s.boundary.gap{kg}.value = v;
    end
end
eq.time_slice{1} = s;

end
