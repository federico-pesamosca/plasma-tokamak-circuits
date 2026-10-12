"""Build a DD 4 machine-description entry for the plasmaless circuit model.

Reads the DD3 DINA em_coupling + pf_active + equilibrium and the ITER_MD
pf_passive used by configure_plasmaless_model.m, converts them to DD 4.1.1 and
writes:
- em_coupling: the three conductor mutual-inductance matrices and the two
  coil-to-grid flux maps (DD3 mutual_grid_active/passive, removed from the DD
  in 3.40.0) as coupling_matrix entries, plus the (R,Z) mesh of the flux maps
  as grid_ggd(1) (one space, nodes only, geometry = [R, Z] of every node),
- pf_active: coil description only (time-dependent data removed),
- pf_passive: loop description,
- wall: limiter outlines used by plot_plasmaless_model.m (ITER_MD 116000/2,
  description_2d(1).limiter.unit(1) = first wall, unit(2) = divertor; only
  name, outline.r and outline.z are kept, in the DD3 point order).
The written entry is reread and compared with the DD3 source.

Mesh of the flux maps (verified against the analytic filament flux of each
coil): the DD3 equilibrium profiles_2d(1) grid of the same DINA entry,
65 uniform R values (grid.dim1) x 129 uniform Z values (grid.dim2). Row k
(0-based) of mutual_grid_* is the node iR, iZ with k = iZ + 129*iR (Z varies
fastest), i.e. mutual_grid_active[:, j].reshape(65, 129) is indexed [iR, iZ].
The grid_ggd nodes are listed in the same order, so coupling_matrix rows_uri(k)
points at object(k). In MATLAB, plasmaless_grid_from_ids(em_coupling) returns
[Mxi, RR, ZZ] (RR = reshape(R_nodes, 129, 65): Z along rows, R along columns;
psi = reshape(Mxi*I, size(RR)) is the flux map, as in plot_plasmaless_model.m);
it reads the entries through get_coupling_matrix, which checks that rows_uri /
columns_uri are the node lists written here (GRID_NODE, pf_active/coil,
pf_passive/loop), the quantity index 1 (magnetic_flux), the grid_ggd(1)
coordinates (r = 4, z = 3) and geometry_content 1 (node_coordinates). Values are full poloidal flux per ampere
(psi = 2*pi*R*A_phi), units H (= Wb/A).

Default output: data/md_dd4 of this repository (the entry shipped with the model,
used by run_plasmaless_imas.m and the MUSCLE3 actor); --output-uri writes elsewhere.
"""

import argparse
from pathlib import Path

import imas
import numpy as np
from imas.util import tree_iter

DD_VERSION = "4.1.1"
DINA_URI = "imas:mdsplus?path=/work/projects/dina/SRO_JINTRAC/15MA_10perc/output"
PF_PASSIVE_URI = "imas:mdsplus?path=/work/imas/shared/imasdb/ITER_MD/3/115005/3"
WALL_URI = "imas:mdsplus?path=/work/imas/shared/imasdb/ITER_MD/3/116000/2"
OUTPUT_URI = f"imas:hdf5?path={Path(__file__).resolve().parents[1] / 'data' / 'md_dd4'}"

# path of the grid nodes (grid_ggd(1), one space, objects of dimension 1 = nodes)
GRID_NODE = "grid_ggd(1)/space(1)/objects_per_dimension(1)/object"
# expected mesh size of the DINA flux maps (nR x nZ)
GRID_SHAPE = (65, 129)

# name, DD3 field, rows_uri node, columns_uri node
# (rows_uri/columns_uri list every index: the DD validation requires their
# size to match the data dimensions, so the implicit "(:)" form is not used)
MATRICES = [
    (
        "mutual_active_active",
        "mutual_active_active",
        "pf_active/coil",
        "pf_active/coil",
    ),
    (
        "mutual_passive_active",
        "mutual_passive_active",
        "pf_passive/loop",
        "pf_active/coil",
    ),
    (
        "mutual_passive_passive",
        "mutual_passive_passive",
        "pf_passive/loop",
        "pf_passive/loop",
    ),
    (
        "mutual_grid_active",
        "mutual_grid_active",
        GRID_NODE,
        "pf_active/coil",
    ),
    (
        "mutual_grid_passive",
        "mutual_grid_passive",
        GRID_NODE,
        "pf_passive/loop",
    ),
]


def read_dd3(uri, name, lazy=False):
    with imas.DBEntry(uri, "r") as entry:
        return entry.get(name, autoconvert=False, lazy=lazy)


def read_grid(uri):
    """(R, Z) node coordinates of the DINA flux maps, in the mutual_grid_* row order.

    The mesh is the equilibrium profiles_2d(1) rectangular grid (dim1 = R,
    dim2 = Z) of the DINA entry; r/z are (nR, nZ) arrays with r[:, 0] == dim1
    and z[0, :] == dim2, so their C-order flattening is k = iZ + nZ*iR.
    """
    with imas.DBEntry(uri, "r") as entry:
        eq = entry.get("equilibrium", autoconvert=False, lazy=True)
        p2d = eq.time_slice[0].profiles_2d[0]
        dim1 = np.array(p2d.grid.dim1)
        dim2 = np.array(p2d.grid.dim2)
        r = np.array(p2d.r)
        z = np.array(p2d.z)
    assert r.shape == z.shape == (dim1.size, dim2.size) == GRID_SHAPE, (
        r.shape,
        z.shape,
    )
    assert np.array_equal(r[:, 0], dim1) and np.array_equal(z[0, :], dim2)
    rr, zz = np.meshgrid(dim1, dim2, indexing="ij")
    assert np.array_equal(rr, r) and np.array_equal(zz, z)
    return dim1, dim2


def set_properties(ids, comment, sources):
    ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_INDEPENDENT
    ids.ids_properties.comment = comment
    ids.ids_properties.provenance.node.resize(1)
    ids.ids_properties.provenance.node[0].path = ""
    ids.ids_properties.provenance.node[0].reference.resize(len(sources))
    for ref, src in zip(ids.ids_properties.provenance.node[0].reference, sources):
        ref.name = src


def set_identifier(node, name, index, description):
    node.name = name
    node.index = index
    node.description = description


def build_grid_ggd(grid_ggd, dim1, dim2, dina_uri):
    """Describe the (R,Z) mesh of the flux maps: one space, nodes only."""
    nr, nz = dim1.size, dim2.size
    # no ggd_identifier option describes a single structured 2-D space of
    # nodes (10 = structured_spaces needs one 1-D space per dimension), so
    # index 0 (unspecified) is used with a free name and a full description
    set_identifier(
        grid_ggd.identifier,
        "DINA_RZ_grid_65x129",
        0,
        f"DINA {nr} x {nz} rectangular (R,Z) grid of the coupling_matrix "
        "mutual_grid_active/mutual_grid_passive flux maps, uniform spacing: "
        f"R = {dim1[0]:.7g}..{dim1[-1]:.7g} m ({nr} values), "
        f"Z = {dim2[0]:.7g}..{dim2[-1]:.7g} m ({nz} values). Nodes are listed "
        f"with Z varying fastest: node k = iZ + {nz}*(iR-1) (1-based k, iR, iZ), "
        f"geometry = [R(iR), Z(iZ)]. Source: equilibrium (DD 3.39.0) "
        f"time_slice(1)/profiles_2d(1)/grid dim1 (R) and dim2 (Z) of {dina_uri}",
    )
    grid_ggd.space.resize(1)
    space = grid_ggd.space[0]
    set_identifier(
        space.identifier,
        "primary_standard",
        1,
        "Primary space defining the standard grid",
    )
    set_identifier(
        space.geometry_type, "standard", 0, "Standard (non-Fourier) geometry"
    )
    space.coordinates_type.resize(2)
    set_identifier(space.coordinates_type[0], "r", 4, "Major radius")
    set_identifier(space.coordinates_type[1], "z", 3, "Vertical coordinate z")
    space.objects_per_dimension.resize(1)
    opd = space.objects_per_dimension[0]
    set_identifier(
        opd.geometry_content, "node_coordinates", 1, "For nodes : node coordinates"
    )
    rr, zz = np.meshgrid(dim1, dim2, indexing="ij")
    nodes = np.column_stack([rr.ravel(), zz.ravel()])  # k = iZ + nz*iR
    opd.object.resize(nodes.shape[0])
    for obj, node in zip(opd.object, nodes):
        obj.geometry = node


def build_em_coupling(em3, dim1, dim2, dina_uri):
    em4 = imas.IDSFactory(DD_VERSION).em_coupling()
    set_properties(
        em4,
        "Conductor mutual inductances and coil-to-grid flux maps (H) copied "
        f"from the DD3 em_coupling mutual_* fields of {dina_uri}; the (R,Z) "
        "mesh of the flux maps (grid_ggd(1)) is the equilibrium profiles_2d "
        "grid of the same entry",
        [f"{dina_uri}#em_coupling", f"{dina_uri}#equilibrium"],
    )
    em4.grid_ggd.resize(1)
    build_grid_ggd(em4.grid_ggd[0], dim1, dim2, dina_uri)
    em4.coupling_matrix.resize(len(MATRICES))
    for cm, (name, field, rows, cols) in zip(em4.coupling_matrix, MATRICES):
        cm.name = name
        cm.quantity.index = 1
        cm.quantity.name = "magnetic_flux"
        cm.quantity.description = "Magnetic flux"
        data = np.array(em3[field])
        if rows == GRID_NODE:
            assert data.shape[0] == dim1.size * dim2.size, (name, data.shape)
        cm.rows_uri = [f"{rows}({i + 1})" for i in range(data.shape[0])]
        cm.columns_uri = [f"{cols}({j + 1})" for j in range(data.shape[1])]
        cm.data = data
    return em4


def build_pf_active(pfa3, dina_uri):
    pfa4 = imas.convert_ids(pfa3, DD_VERSION)
    # remove all time-dependent data (coil/supply waveforms, code output, time)
    for node in tree_iter(pfa4, leaf_only=True):
        if node.metadata.type.is_dynamic and node.has_value:
            if node.metadata.ndim == 0:
                raise ValueError(f"Cannot clear dynamic scalar {node.metadata.path}")
            node.value = np.array([], dtype=node.value.dtype)
    set_properties(
        pfa4,
        f"Coil description of {dina_uri} (DD3 pf_active converted to DD "
        f"{DD_VERSION}, time-dependent data removed)",
        [f"{dina_uri}#pf_active"],
    )
    return pfa4


def build_pf_passive(pfp3, pf_passive_uri):
    pfp4 = imas.convert_ids(pfp3, DD_VERSION)
    set_properties(
        pfp4,
        f"Passive loops of {pf_passive_uri} converted to DD {DD_VERSION}",
        [f"{pf_passive_uri}#pf_passive"],
    )
    return pfp4


def build_wall(wall3, wall_uri):
    """Wall IDS reduced to the two limiter units (first wall, divertor) outlines."""
    wall4 = imas.IDSFactory(DD_VERSION).wall()
    set_properties(
        wall4,
        f"Limiter units 1 (first wall) and 2 (divertor) outlines of {wall_uri} "
        f"(DD3 wall converted to DD {DD_VERSION}; other content not copied)",
        [f"{wall_uri}#wall"],
    )
    wall4.description_2d.resize(1)
    units3 = wall3.description_2d[0].limiter.unit
    assert len(units3) >= 2, len(units3)
    units4 = wall4.description_2d[0].limiter.unit
    units4.resize(2)
    for u3, u4 in zip(list(units3)[:2], units4):
        u4.name = u3.name.value
        u4.outline.r = np.array(u3.outline.r)
        u4.outline.z = np.array(u3.outline.z)
    return wall4


def check_wall(wall4, wall3):
    units3 = wall3.description_2d[0].limiter.unit
    units4 = wall4.description_2d[0].limiter.unit
    assert len(units4) == 2
    for k in range(2):
        for f in ("r", "z"):
            assert np.array_equal(
                np.asarray(getattr(units4[k].outline, f)),
                np.asarray(getattr(units3[k].outline, f)),
            ), (k, f)
    print(
        "  wall limiter units: "
        f"{[(u.name.value, len(u.outline.r)) for u in units4]}"
    )


def check(output_uri, em3, pfa3, pfp3, dim1, dim2, wall3):
    with imas.DBEntry(output_uri, "r") as entry:
        em4 = entry.get("em_coupling")
        pfa4 = entry.get("pf_active")
        pfp4 = entry.get("pf_passive")
        wall4 = entry.get("wall")
    check_wall(wall4, wall3)
    names = [cm.name.value for cm in em4.coupling_matrix]
    for name, field, rows, cols in MATRICES:
        cm = em4.coupling_matrix[names.index(name)]
        data = np.asarray(cm.data)
        assert np.array_equal(data, np.asarray(em3[field])), name
        assert list(cm.rows_uri) == [
            f"{rows}({i + 1})" for i in range(data.shape[0])
        ], name
        assert list(cm.columns_uri) == [
            f"{cols}({j + 1})" for j in range(data.shape[1])
        ], name
    # grid_ggd nodes reproduce the equilibrium grid in the mutual_grid_* row order
    objects = em4.grid_ggd[0].space[0].objects_per_dimension[0].object
    nodes = np.array([np.asarray(obj.geometry) for obj in objects])
    assert nodes.shape == (dim1.size * dim2.size, 2), nodes.shape
    r_nodes = nodes[:, 0].reshape(dim1.size, dim2.size)  # C order: [iR, iZ]
    z_nodes = nodes[:, 1].reshape(dim1.size, dim2.size)
    assert np.array_equal(r_nodes, np.repeat(dim1[:, None], dim2.size, axis=1))
    assert np.array_equal(z_nodes, np.repeat(dim2[None, :], dim1.size, axis=0))
    assert len(pfa4.coil) == len(pfa3.coil)
    for c3, c4 in zip(pfa3.coil, pfa4.coil):
        assert c4.name.value == c3.identifier.value, (c4.name, c3.identifier)
        assert c4.resistance.value == c3.resistance.value, c4.name
        assert len(c4.current.data) == 0 and len(c4.voltage.data) == 0
    assert len(pfp4.loop) == len(pfp3.loop)
    for l3, l4 in zip(pfp3.loop, pfp4.loop):
        assert l4.name.value == l3.name.value
        assert l4.resistance.value == l3.resistance.value, l4.name
    print(f"Checked {output_uri} (DD {em4.ids_properties.version_put.data_dictionary})")
    for cm in em4.coupling_matrix:
        print(f"  em_coupling {cm.name.value}: {cm.data.shape}")
    print(
        f"  em_coupling grid_ggd(1): {nodes.shape[0]} nodes = {dim1.size} R x {dim2.size} Z, "
        f"R = {dim1[0]:.7g}..{dim1[-1]:.7g} m, Z = {dim2[0]:.7g}..{dim2[-1]:.7g} m "
        "(node k = iZ + nZ*iR, Z fastest)"
    )
    print(f"  pf_active coils: {[c.name.value for c in pfa4.coil]}")
    print(f"  pf_active resistances: {[float(c.resistance) for c in pfa4.coil]}")
    print(f"  pf_passive loops: {len(pfp4.loop)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dina-uri", default=DINA_URI)
    parser.add_argument("--pf-passive-uri", default=PF_PASSIVE_URI)
    parser.add_argument("--wall-uri", default=WALL_URI)
    parser.add_argument("--output-uri", default=OUTPUT_URI)
    args = parser.parse_args()

    em3 = read_dd3(args.dina_uri, "em_coupling")
    pfa3 = read_dd3(args.dina_uri, "pf_active")
    pfp3 = read_dd3(args.pf_passive_uri, "pf_passive")
    wall3 = read_dd3(args.wall_uri, "wall")
    dim1, dim2 = read_grid(args.dina_uri)

    em4 = build_em_coupling(em3, dim1, dim2, args.dina_uri)
    pfa4 = build_pf_active(pfa3, args.dina_uri)
    pfp4 = build_pf_passive(pfp3, args.pf_passive_uri)
    wall4 = build_wall(wall3, args.wall_uri)

    with imas.DBEntry(args.output_uri, "w", dd_version=DD_VERSION) as entry:
        entry.put(em4)
        entry.put(pfa4)
        entry.put(pfp4)
        entry.put(wall4)

    check(args.output_uri, em3, pfa3, pfp3, dim1, dim2, wall3)


if __name__ == "__main__":
    main()
