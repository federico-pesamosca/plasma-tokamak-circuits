function [Mxi, RR, ZZ, Mxa, Mxe] = plasmaless_grid_from_ids(em_coupling)
% coil-to-grid flux maps and their (R,Z) mesh from a DD4 em_coupling
% em_coupling: DD4 em_coupling with the coupling_matrix entries
%              mutual_grid_active (Ngrid x n_coils) and mutual_grid_passive
%              (Ngrid x n_loops) and the mesh in grid_ggd{1} (one space, nodes
%              only), as written by tools/convert_md_dd3_to_dd4.py
% Mxi:  [Mxa, Mxe] (Ngrid x (n_coils+n_loops)), full poloidal flux per ampere
%       (psi = 2*pi*R*A_phi) at every grid node, units H (= Wb/A); columns in
%       the pf_active coil order then the pf_passive loop order
% RR, ZZ: nZ x nR node coordinates (m), Z along rows, R along columns, so that
%       psi = reshape(Mxi*[Ia; Ie], size(RR)) is the flux map on the mesh and
%       contour(RR, ZZ, psi) plots it, as in plot_plasmaless_model.m
% Mxa, Mxe: the two maps separately
%
% Node ordering: the grid_ggd nodes (and the rows of the maps) are listed with
% Z varying fastest, node k = iZ + nZ*(iR-1) (1-based), geometry = [R(iR); Z(iZ)];
% reshape(..., nZ, nR) therefore puts Z along rows (MATLAB is column-major).
% The mesh must be a tensor-product grid (every column of RR constant, every
% row of ZZ constant); otherwise an error is raised.
%
% DD3 entries (DINA, DD <= 3.39): the maps are the fields mutual_grid_active and
% mutual_grid_passive of em_coupling (removed from the DD in 3.40.0, so they need
% an IMAS module with DD <= 3.39), and the mesh is the equilibrium profiles_2d{1}
% grid of the same entry (RR = r', ZZ = z'); this function does not read them.

node = 'grid_ggd(1)/space(1)/objects_per_dimension(1)/object';
Mxa = get_coupling_matrix(em_coupling, 'mutual_grid_active', node, 'pf_active/coil');
Mxe = get_coupling_matrix(em_coupling, 'mutual_grid_passive', node, 'pf_passive/loop');
if size(Mxa,1) ~= size(Mxe,1)
    error('plasmaless:grid','mutual_grid_active has %d rows, mutual_grid_passive %d.', ...
        size(Mxa,1), size(Mxe,1));
end
Mxi = [Mxa, Mxe];
n = size(Mxi,1);

% mesh: grid_ggd{1}, one space with coordinates (r, z), objects of dimension 1 = nodes
if ~isfield(em_coupling,'grid_ggd') || numel(em_coupling.grid_ggd) < 1
    error('plasmaless:grid','em_coupling has no grid_ggd (mesh of the flux maps).');
end
grid = em_coupling.grid_ggd{1};
if numel(grid.space) ~= 1
    error('plasmaless:grid','grid_ggd{1} has %d spaces, expected 1.', numel(grid.space));
end
space = grid.space{1};
ct = space.coordinates_type;
if numel(ct) ~= 2 || ct{1}.index ~= 4 || ct{2}.index ~= 3
    error('plasmaless:grid', ...
        'grid_ggd{1}.space{1} coordinates are not (r, z) (coordinate_identifier 4, 3).');
end
if numel(space.objects_per_dimension) < 1
    error('plasmaless:grid','grid_ggd{1}.space{1} has no objects_per_dimension.');
end
opd = space.objects_per_dimension{1};
if opd.geometry_content.index ~= 1
    error('plasmaless:grid', ...
        'grid_ggd{1} nodes: geometry_content index %d is not 1 (node_coordinates).', ...
        opd.geometry_content.index);
end
if numel(opd.object) ~= n
    error('plasmaless:grid','grid_ggd{1} has %d nodes, the flux maps %d rows.', ...
        numel(opd.object), n);
end
R = zeros(n,1);
Z = zeros(n,1);
for k = 1:n
    geom = opd.object{k}.geometry;
    if numel(geom) ~= 2
        error('plasmaless:grid','grid_ggd{1} node %d: geometry has %d values, expected [R; Z].', ...
            k, numel(geom));
    end
    R(k) = geom(1);
    Z(k) = geom(2);
end

% tensor-product grid with Z varying fastest: nZ = length of the first run of constant R
nZ = find(R ~= R(1), 1) - 1;
if isempty(nZ) || mod(n, nZ) ~= 0
    error('plasmaless:grid','grid_ggd{1} nodes are not a rectangular (R,Z) grid with Z varying fastest.');
end
nR = n / nZ;
RR = reshape(R, nZ, nR);
ZZ = reshape(Z, nZ, nR);
if any(RR ~= RR(1,:), 'all') || any(ZZ ~= ZZ(:,1), 'all')
    error('plasmaless:grid', ...
        'grid_ggd{1} nodes are not a tensor-product (R,Z) grid (R must depend on the column only, Z on the row only).');
end

end
