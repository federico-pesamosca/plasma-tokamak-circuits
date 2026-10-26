"""Voltage driver for the standalone MUSCLE3 test of the plasmaless actor.

Stands in for the PDS controller (plasmaless_standalone.ymmsl): sends the F_INIT
equilibrium and pf_active, then answers every actor output with the power-supply
voltages of the pulse train of timeseries_plasmaless_model.m (one +/-10 V pulse per
supply, in the Vps order CS3U, CS2U, CS1U, CS1L, CS2L, CS3L, PF1 .. PF5, VS1, PF6, VS3),
as run_plasmaless_imas.m applies them (cold start, u(t_k) over [t_k, t_k+dt]). Stores
the received currents and the observables carried by equilibrium_in (Y, NaN where
absent) and saves them in the run_plasmaless_imas.m results layout (time, Ia, Ie, Y,
coil_names, loop_names, y_names), plus Vps (the voltages this driver sent, nt x 14,
one row per time; the last row is the pulse-train value at t_end, not sent) and
input_names (the 14 supply keys in Vps order).

Setting `command` (default 'supply') selects the format of the voltage command, the
two inputs the actor (and NICE evo) accept:
  'supply': the 14 power-supply voltages in pf_active.supply(:).voltage.data, no
            names and no coils (as a controller would send them; the actor names the
            supplies by index). The F_INIT pf_active carries the coils (currents 0,
            voltages 0, the actor's contract) and the 14 supplies copied from the
            machine description with the pulse-train value at t0. This mode applies
            the whole pulse train, VS1 included: the results equal those of
            run_plasmaless_imas.m with the same dt.
  'coil':   the 14 coil voltages in pf_active.coil(:).voltage.data, no names (the PDS
            controller's format, as before): CS/PF supplies one-to-one (CS1U <-
            Vps_CS1U, CS1L <- Vps_CS1L), VS1 dropped, VS3 split as VS3L = +Vps_VS3,
            VS3U = 0 so that the actor's fallback V_VS3L - V_VS3U = Vps_VS3. The
            F_INIT pf_active has no supplies and the same coil-sense voltages at t0.
            This mode loses the VS1 pulse (t = 22 .. 23 s of the pulse train, the
            12th supply): the PF2 .. PF5 currents differ from run_plasmaless_imas.m
            there; everything else is identical.

Ports: O_I equilibrium_out, pf_active_out (to the actor's F_INIT), voltages_out (to
pf_active_s); S equilibrium_in, pf_active_in, pf_passive_in (actor's O_I).
The F_INIT data goes out before the O_I/S loop, without voltages_out: libmuscle
warns once that the first receive does not adhere to the MMSF (harmless here).
Settings: dt, t_end (shared with the actor), md_uri, results_file (absolute paths,
set by run_plasmaless_muscle3.sh in a settings file of the run directory), command
(optional, 'supply' | 'coil', see above; run_plasmaless_muscle3.sh writes it from
PLASMALESS_M3_COMMAND).
"""

import imas
import numpy as np
import scipy.io
from libmuscle import Instance, Message
from ymmsl import Operator

TSTART = 0.0  # run_plasmaless_imas.m:17
# observables vocabulary of imas_model/plasmaless_observables_to_equilibrium.m
Y_NAMES = ["rmag", "zmag", "rX", "zX", "Ip", "gap1", "gap2", "gap3", "Dist2ndSep", "rc", "zc"]
# the 14 power supplies in the order of the pulse train (Vps order of
# imas_model/plasmaless_circuit_topology.m, = the machine-description supply order)
INPUT_NAMES = ["CS3U", "CS2U", "CS1U", "CS1L", "CS2L", "CS3L", "PF1", "PF2", "PF3",
               "PF4", "PF5", "VS1", "PF6", "VS3"]


def _value(node):
    """Float of a filled IDS scalar, NaN when it is not filled."""
    return float(node) if node.has_value else np.nan


def supply_key(name):
    """Short supply key of a machine-description supply name ('CS3U-1-2' -> 'CS3U',
    'VS1-1-2-3-4-5-6' -> 'VS1'): the name up to the first '-', as
    imas_model/plasmaless_supply_key.m."""
    return str(name).split("-", 1)[0].strip()


def observables(ts):
    """Row of Y_NAMES values of one equilibrium time slice (DD4 fields written by
    plasmaless_observables_to_equilibrium.m): magnetic axis, the first contour_tree
    node with critical_type 1 (X-point), ip, up to 3 gaps, current centre r/z
    (vertical model); Dist2ndSep has no DD4 field and is always NaN."""
    row = np.full(len(Y_NAMES), np.nan)
    row[0] = _value(ts.global_quantities.magnetic_axis.r)
    row[1] = _value(ts.global_quantities.magnetic_axis.z)
    for node in ts.contour_tree.node:
        if node.critical_type.has_value and int(node.critical_type) == 1:
            row[2] = _value(node.r)
            row[3] = _value(node.z)
            break
    row[4] = _value(ts.global_quantities.ip)
    for k, gap in enumerate(ts.boundary.gap):
        if k >= 3:
            break
        row[5 + k] = _value(gap.value)
    row[9] = _value(ts.global_quantities.current_centre.r)
    row[10] = _value(ts.global_quantities.current_centre.z)
    return row


def pulse_train(t):
    """Power-supply voltages Vps at the times t (timeseries_plasmaless_model.m:4-29)."""
    n, a, pulse_width, break_time = 14, 10.0, 1.0, 1.0  # .m:4-7
    u = np.zeros((len(t), n))
    for k in range(1, n + 1):  # .m:18-29
        t_start = (k - 1) * (pulse_width + break_time)
        t_mid = t_start + pulse_width / 2
        t_end = t_start + pulse_width
        u[(t >= t_start) & (t < t_mid), k - 1] = a
        u[(t >= t_mid) & (t < t_end), k - 1] = -a
    return u


def voltages_at(t, dt):
    """Pulse train at the grid point TSTART + k*dt closest to t.

    The actor accumulates t += dt (round-off up to ~1e-12 s), which flips pulse edges
    t >= t_start at a few steps; the grid point gives the same voltages as the
    (tstart:dt:tend) samples of the .m file (checked for dt = 0.01 and 0.001).
    """
    k = round((t - TSTART) / dt)
    return pulse_train(np.array([TSTART + k * dt]))[0]


def coil_voltages(u, coil_names):
    """Coil-sense image of the supply voltages u (Vps order), 'coil' command mode:
    CS/PF coils one-to-one with their supply (CS1U <- Vps_CS1U, CS1L <- Vps_CS1L),
    VS1 dropped (no coil carries it), VS3L = +Vps_VS3 and VS3U = 0 (the actor's
    fallback takes the VS3 supply voltage as V_VS3L - V_VS3U)."""
    v = np.zeros(len(coil_names))
    for i, name in enumerate(coil_names):
        if name == "VS3L":
            v[i] = u[INPUT_NAMES.index("VS3")]
        elif name == "VS3U":
            v[i] = 0.0
        else:
            v[i] = u[INPUT_NAMES.index(name)]
    return v


inst = Instance(
    {
        Operator.O_I: ["equilibrium_out", "pf_active_out", "voltages_out"],
        Operator.S: ["equilibrium_in", "pf_active_in", "pf_passive_in"],
    }
)
while inst.reuse_instance():
    dt = inst.get_setting("dt", "float")
    t_end = inst.get_setting("t_end", "float")
    md_uri = inst.get_setting("md_uri", "str")
    results_file = inst.get_setting("results_file", "str")
    # libmuscle 0.10: default= (keyword-only, since 0.9.0) is returned when the
    # setting is not set; without it a missing setting raises KeyError
    command = inst.get_setting("command", "str", default="supply")
    if command not in ("supply", "coil"):
        raise ValueError(f"voltage_driver: setting command must be 'supply' or 'coil' (got '{command}')")

    with imas.DBEntry(md_uri, "r") as entry:
        pfa_md = entry.get("pf_active")
        factory = entry.factory

    # the machine-description supplies must be the 14 of the pulse train, in its order
    md_keys = [supply_key(s.name) for s in pfa_md.supply]
    if md_keys != INPUT_NAMES:
        raise ValueError(
            f"voltage_driver: machine-description supplies {md_keys} are not the pulse-train "
            f"supplies {INPUT_NAMES} (names or order)"
        )
    md_coil_names = [supply_key(c.name) for c in pfa_md.coil]

    # F_INIT equilibrium: only gives t0 (message timestamp) and the reference that
    # the actor passes through. There is no plasma: ip and geometric axis are 0.
    # It ends at t_end + dt because the actor's accumulated time can exceed t_end
    # by round-off and the actor interpolates the reference at that time.
    eq = factory.equilibrium()
    eq.ids_properties.homogeneous_time = 1
    eq.time = np.array([TSTART, t_end + dt])
    eq.time_slice.resize(2)
    for ts, t in zip(eq.time_slice, eq.time):
        ts.time = t
        ts.global_quantities.ip = 0.0
        ts.boundary.geometric_axis.r = 0.0
        ts.boundary.geometric_axis.z = 0.0
    inst.send("equilibrium_out", Message(TSTART, None, eq.serialize()))

    # F_INIT pf_active: coils of the machine description (same names and order as
    # the model), cold start (currents 0) and the pulse-train voltages at t0:
    # 'supply': coil voltages 0 (contract) and the 14 supplies of the machine
    # description with u0; 'coil': no supplies, coil voltages = coil image of u0
    u0 = voltages_at(TSTART, dt)
    pfa = factory.pf_active()
    pfa.ids_properties.homogeneous_time = 1
    pfa.time = np.array([TSTART])
    pfa.coil.resize(len(pfa_md.coil))
    v0 = coil_voltages(u0, md_coil_names) if command == "coil" else np.zeros(len(pfa_md.coil))
    for i, (c, c_md) in enumerate(zip(pfa.coil, pfa_md.coil)):
        c.name = c_md.name
        c.current.data = np.array([0.0])
        c.voltage.data = np.array([v0[i]])
    if command == "supply":
        pfa.supply.resize(len(pfa_md.supply))
        for j, (s, s_md) in enumerate(zip(pfa.supply, pfa_md.supply)):
            s.name = s_md.name
            s.voltage.data = np.array([u0[j]])
    inst.send("pf_active_out", Message(TSTART, None, pfa.serialize()))

    time, ia, ie, vps = [TSTART], [np.zeros(len(pfa_md.coil))], None, [u0]
    y = [observables(eq.time_slice[0])]  # F_INIT values at t0 (ip 0, rest NaN)
    while True:
        msg_eq = inst.receive("equilibrium_in")
        msg_pfa = inst.receive("pf_active_in")
        msg_pfp = inst.receive("pf_passive_in")
        eq_out = factory.equilibrium()
        eq_out.deserialize(msg_eq.data)
        y.append(observables(eq_out.time_slice[0]))
        pfa_out = factory.pf_active()
        pfa_out.deserialize(msg_pfa.data)
        pfp_out = factory.pf_passive()
        pfp_out.deserialize(msg_pfp.data)
        if ie is None:  # vessel currents 0 at t0 (actor assumption)
            ie = [np.zeros(len(pfp_out.loop))]
            coil_names = [str(c.name) for c in pfa_out.coil]
            loop_names = [str(lp.name) for lp in pfp_out.loop]
        t = msg_pfa.timestamp
        time.append(t)
        ia.append([float(c.current.data[0]) for c in pfa_out.coil])
        ie.append([float(lp.current[0]) for lp in pfp_out.loop])
        # voltages applied over [t, t+dt] (pulse-train value at t; the last one,
        # at t_end, is recorded but not sent)
        u = voltages_at(t, dt)
        vps.append(u)
        if msg_pfa.next_timestamp is None:
            break
        # no names (as the PDS controller; the actor sets them by index)
        cmd = factory.pf_active()
        cmd.ids_properties.homogeneous_time = 1
        cmd.time = np.array([t])
        if command == "supply":
            cmd.supply.resize(len(u))
            for s, v in zip(cmd.supply, u):
                s.voltage.data = np.array([v])
        else:
            v_coil = coil_voltages(u, md_coil_names)
            cmd.coil.resize(len(v_coil))
            for c, v in zip(cmd.coil, v_coil):
                c.voltage.data = np.array([v])
        inst.send("voltages_out", Message(t, None, cmd.serialize()))

    # same variables and shapes as run_plasmaless_imas.m:102, plus Vps / input_names
    scipy.io.savemat(
        results_file,
        {
            "time": np.array(time),
            "Ia": np.array(ia),
            "Ie": np.array(ie),
            "Y": np.array(y),
            "Vps": np.array(vps),
            "coil_names": np.array(coil_names, dtype=object).reshape(1, -1),
            "loop_names": np.array(loop_names, dtype=object).reshape(1, -1),
            "y_names": np.array(Y_NAMES, dtype=object).reshape(1, -1),
            "input_names": np.array(INPUT_NAMES, dtype=object).reshape(1, -1),
        },
        oned_as="column",
    )
    print(
        f"voltage_driver: command '{command}', {len(time) - 1} steps, t = {time[0]} .. {time[-1]}, "
        f"results saved to {results_file}"
    )
