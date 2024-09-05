import numpy as np
from sympy import symbols
from sympy.physics.mechanics import (
    Particle, Point,
    ReferenceFrame, RigidBody,
    System, dynamicsymbols
)
from pathway import LinearPathway


def particle_pend(n, force=None):
    """n defines the number of small masses we want to add in the wire."""
    L = symbols('L', positive=True)

    # Creating Reference frames, angles, coordinates, and speeds
    frames = symbols(f'N1:{n + 2}', cls=ReferenceFrame)
    angles = dynamicsymbols(f'alpha1:{n + 2}')
    betas = symbols(f'beta1:{n + 2}')
    qs = dynamicsymbols(f'q1:{2 * (n + 1) + 1}')
    us = dynamicsymbols(f'u1:{2 * (n + 1) + 1}')
    qc = dynamicsymbols('qc1:4')
    masses = symbols(f'm1:{n + 1}')
    M = symbols('M')
    masses = list(masses) + [M]

    metric = [1, 1, 0]

    t, g = symbols('t g')
    l1 = symbols('l1', positive=True)
    natural_lengths = symbols(f'l0_1:{n + 2}')
    ks = symbols(f'k1:{n + 2}')

    # Creating the Points
    ps = symbols(f'P1:{n + 2}', cls=Point)

    system = System()

    # Setting these points in space
    j = 0
    for i in range(n):
        ps[i].set_pos(
            system.fixed_point,
            qs[j] * metric[0] * system.frame.x +
            qs[j + 1] * metric[1] * system.frame.y
        )
        j += 2

    ps[-1].set_pos(
        system.fixed_point,
        (qs[-2] + qc[0]) * metric[0] * system.frame.x +
        (qs[-1] + qc[1]) * metric[1] * system.frame.y
    )

    # Rotating the frames with respect to system frame
    for i, j in zip(frames, angles):
        i.orient_axis(system.frame, j, system.frame.z)

    # Creating Particles on the wire
    bodies = [
        Particle(f'B{h}', mass=i, point=j)
        for h, i, j in zip(range(1, n + 1), masses, ps[:-1])
    ]

    # Creating the RigidBody at the bottom of the wire
    last_body = RigidBody(
        f'B{n + 1}', mass=masses[-1], masscenter=ps[-1], frame=frames[-1]
    )

    bodies = bodies + [last_body]

    # Adding bodies to the system
    for body in bodies:
        system.add_bodies(body)

    # Assigning velocities to the Particles
    j = 0
    for i in range(n + 1):
        bodies[i].masscenter.set_vel(
            system.frame,
            us[j] * system.frame.x + us[j + 1] * system.frame.y
        )
        j += 2

    # Assigning velocities to the RigidBody
    bodies[-1].masscenter.set_vel(
        system.frame, us[-2] * system.frame.x + us[-1] * system.frame.y
    )
    bodies[-1].frame.set_ang_vel(frames[-1], 0)

    # Adding coordinates, speeds, and kdes to the system
    for k in zip(qs, us):
        system.add_coordinates(k[0])
        system.add_speeds(k[1])
        system.add_kdes(k[1] - k[0].diff(t))

    # Applying gravity
    system.apply_uniform_gravity(-g * system.frame.y)

    # Creating linear pathways between the points and applying force along them
    linear_points = [system.fixed_point] + list(ps)

    linear_paths = []
    for i in range(1, len(linear_points) - 1):
        lin_paths_a = LinearPathway(linear_points[i - 1], linear_points[i])
        lin_paths_b = LinearPathway(linear_points[i + 1], linear_points[i])
        linear_paths.append(lin_paths_a)
        linear_paths.append(lin_paths_b)
    linear_paths.append(LinearPathway(linear_points[-2], linear_points[-1]))

    deltas_ = []
    delta_iterator = 0
    for paths in linear_paths[::2]:
        delta_val = (
            paths.attachments[0].pos_from(paths.attachments[1]).magnitude() -
            natural_lengths[delta_iterator]
        )
        deltas_.append(delta_val)
        delta_iterator += 1

    forces_values = [i * j for i, j in zip(ks, deltas_)]

    forces_a = [forces_values[0]] + [
        val for val in forces_values[1:] for _ in (0, 1)
    ]
    force_iterator = 0
    for path in linear_paths:
        system.add_loads(path.to_loads(-forces_a[force_iterator])[1])
        force_iterator += 1

    # Validating the system
    system.validate_system()

    # Getting the Kane's equations
    system.form_eoms()
    kanel = system.eom_method.to_linearizer()

    # The operating point is x coordinates are 0
    # y coordinates are length of wire so l1, 2l1, 3l1 and so on
    # the additional coordinates of last body go to zero
    op_vals = (
        {qs[i]: 0 for i in range(0, 2 * (n + 1), 2)} |
        {qs[i]: (i - i // 2) * l1 for i in range(1, 2 * (n + 1), 2)} |
        {qc[i]: 0 for i in range(0, 3, 1)}
    )

    A, B = kanel.linearize(A_and_B=True, op_point=op_vals)

    # Substituting the value of natural lengths in terms of l1, masses,
    # and spring constants
    lsdict = {
        i: l1 - (np.sum(masses[m:]) * g) / k
        for i, m, k in zip(natural_lengths, range(n + 2), ks)
    }

    A = A.subs(lsdict)
    B = B.subs(lsdict)

    A = A.subs({l1: L / (n + 1)})

    return A, B, system
