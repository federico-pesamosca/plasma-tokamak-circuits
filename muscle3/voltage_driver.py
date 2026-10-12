"""Voltage driver for the standalone MUSCLE3 test of the plasmaless actor.

Stands in for the PDS controller (plasmaless_standalone.ymmsl): sends the F_INIT
equilibrium and pf_active, then answers every actor output with the coil voltages of
the pulse train of timeseries_plasmaless_model.m, as run_plasmaless_imas.m applies
them (cold start, u(t_k) over [t_k, t_k+dt]). Stores the received currents and saves
them in the run_plasmaless_imas.m results layout (time, Ia, Ie, coil_names,
loop_names), same variables as run_plasmaless_imas.m.

Ports: O_I equilibrium_out, pf_active_out (to the actor's F_INIT), voltages_out (to
pf_active_s); S equilibrium_in, pf_active_in, pf_passive_in (actor's O_I).
The F_INIT data goes out before the O_I/S loop, without voltages_out: libmuscle
warns once that the first receive does not adhere to the MMSF (harmless here).
Settings: dt, t_end (shared with the actor), md_uri, results_file (absolute paths,
set by run_plasmaless_muscle3.sh in a settings file of the run directory).
"""

import imas
import numpy as np
import scipy.io
from libmuscle import Instance, Message
from ymmsl import Operator

TSTART = 0.0  # run_plasmaless_imas.m:17


def pulse_train(t):
    """Coil voltages at the times t (timeseries_plasmaless_model.m:4-29)."""
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

    with imas.DBEntry(md_uri, "r") as entry:
        pfa_md = entry.get("pf_active")
        factory = entry.factory

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
    # the model), cold start (currents 0) and the pulse-train voltages at t0
    u0 = voltages_at(TSTART, dt)
    pfa = factory.pf_active()
    pfa.ids_properties.homogeneous_time = 1
    pfa.time = np.array([TSTART])
    pfa.coil.resize(len(pfa_md.coil))
    for i, (c, c_md) in enumerate(zip(pfa.coil, pfa_md.coil)):
        c.name = c_md.name
        c.current.data = np.array([0.0])
        c.voltage.data = np.array([u0[i]])
    inst.send("pf_active_out", Message(TSTART, None, pfa.serialize()))

    time, ia, ie = [TSTART], [np.zeros(len(pfa_md.coil))], None
    while True:
        msg_eq = inst.receive("equilibrium_in")
        msg_pfa = inst.receive("pf_active_in")
        msg_pfp = inst.receive("pf_passive_in")
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
        if msg_pfa.next_timestamp is None:
            break
        # voltages applied over [t, t+dt]; no coil names (as the PDS controller)
        u = voltages_at(t, dt)
        cmd = factory.pf_active()
        cmd.ids_properties.homogeneous_time = 1
        cmd.time = np.array([t])
        cmd.coil.resize(len(u))
        for c, v in zip(cmd.coil, u):
            c.voltage.data = np.array([v])
        inst.send("voltages_out", Message(t, None, cmd.serialize()))

    # same variables and shapes as run_plasmaless_imas.m:102
    scipy.io.savemat(
        results_file,
        {
            "time": np.array(time),
            "Ia": np.array(ia),
            "Ie": np.array(ie),
            "coil_names": np.array(coil_names, dtype=object).reshape(1, -1),
            "loop_names": np.array(loop_names, dtype=object).reshape(1, -1),
        },
        oned_as="column",
    )
    print(
        f"voltage_driver: {len(time) - 1} steps, t = {time[0]} .. {time[-1]}, "
        f"results saved to {results_file}"
    )
