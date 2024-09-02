from sympy import zeros, symbols
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
from sympy.physics.mechanics import Point, RigidBody, System, Force
from pathway import LinearPathway
import numpy as np


def make_wire_suspy(n, linpath, index, global_frame, suspension_body):
    """
    Function to create a wire system with a specified number of small masses
    suspended between two points.

    Parameters:
    n (int): Number of small masses to add to the wire.
    linpath (dict): Dictionary containing path information.
    index (int): Index of the pathway.
    global_frame (ReferenceFrame): The global reference frame.
    suspension_body (RigidBody): The suspension body.

    Returns:
    tuple: A tuple containing the system and matrix A.
    """
    point1 = linpath['paths'][index].attachments[0]
    point2 = linpath['paths'][index].attachments[1]

    bodies_str = ''.join([point1.name[0], point2.name[0]])

    body1 = linpath['bodies'][0]
    body2 = linpath['bodies'][-1]

    mass1 = linpath['bodies'][0].mass
    mass2 = linpath['bodies'][1].mass

    gb_frame = global_frame
    ind_str = index + 1
    l1 = symbols('l1', positive=True)
    L = symbols(f'L_{bodies_str}{ind_str}{ind_str}', positive=True)
    F, g, t = symbols('F, g, t')

    frames = symbols(f'N{ind_str}{ind_str}_{bodies_str}1:{n+2}', cls=ReferenceFrame)
    angles = dynamicsymbols(f'alpha{ind_str}{ind_str}_{bodies_str}1:{n+2}')
    betas = symbols(f'beta{ind_str}{ind_str}_{bodies_str}1:{n+2}')
    qs = dynamicsymbols(f'q{ind_str}{ind_str}_{bodies_str}1:{3*(n+1)+1}')
    us = dynamicsymbols(f'u{ind_str}{ind_str}_{bodies_str}1:{3*(n+1)+1}')

    masses = symbols(f'm{ind_str}{ind_str}_{bodies_str}1:{n+2}')
    M = symbols('M')

    model = suspension_body.model

    op_vals = (
        {qs[i]: 0 for i in range(0, 3*(n+1), 3)}
        | {qs[i]: 0 for i in range(1, 3*(n+1), 3)}
        | {qs[i]: -(i - i//3)*l1/2 for i in range(2, 3*(n+1), 3)}
        | {F: mass1 * g}
    )

    dist21 = point2.pos_from(point1).express(gb_frame)
    distx = dist21.dot(gb_frame.x)
    disty = dist21.dot(gb_frame.y)
    distz = dist21.dot(gb_frame.z)

    angle = distz / dist21.magnitude()
    angles[-1] = angle - np.sum([i for i in angles[:-1]])

    S = suspension_body
    wall = RigidBody('W', masscenter=S.global_origin, frame=gb_frame)
    system = System.from_newtonian(wall)

    for i, j in zip(frames, angles):
        i.orient_axis(system.frame, j, system.frame.z)

    ps = symbols(f'P{ind_str}{ind_str}_{bodies_str}1:{n+2}', cls=Point)

    j = 0
    for i in range(n+1):
        ps[i].set_pos(point1, qs[j]*gb_frame.x + qs[j+1]*gb_frame.y + qs[j+2]*gb_frame.z)
        j += 3

    op_point_wire = {
        qs[-3]: distx, qs[-2]: disty, qs[-1]: distz
    }

    bodies = [
        RigidBody(
            f'B{ind_str}{ind_str}_{bodies_str}{h}', mass=i*j,
            masscenter=k, frame=l
        )
        for h, i, j, k, l in zip(range(1, n+2), betas, masses, ps, frames)
    ]

    for body in bodies:
        system.add_bodies(body)

    j = 0
    for i in range(n+1):
        bodies[i].masscenter.set_vel(
            system.frame,
            us[j]*gb_frame.x + us[j+1]*gb_frame.y + us[j+2]*gb_frame.z
        )
        bodies[i].frame.set_ang_vel(frames[i], 0)
        j += 3

    for q, u in zip(qs, us):
        system.add_coordinates(q)
        system.add_speeds(u)
        system.add_kdes(u - q.diff(t))

    system.apply_uniform_gravity(-g * wall.z)

    linear_points = [point1] + [i for i in ps]

    linear_paths = []
    for i in range(1, len(linear_points) - 1):
        lin_paths_a = LinearPathway(linear_points[i - 1], linear_points[i])
        lin_paths_b = LinearPathway(linear_points[i + 1], linear_points[i])
        linear_paths.append(lin_paths_a)
        linear_paths.append(lin_paths_b)

    linear_paths.append(LinearPathway(linear_points[-2], linear_points[-1]))

    for path in linear_paths[:-1]:
        system.add_loads(path.to_loads(-F)[1])

    system.validate_system()

    system.form_eoms()
    kanel = system.eom_method.to_linearizer()

    A, B = kanel.linearize(A_and_B=True, op_point=op_vals)

    A = A.subs({l1: L/(n+1)})
    B = B.subs({l1: L/(n+1)})

    return system, A

