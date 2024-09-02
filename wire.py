from sympy import symbols
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point, RigidBody, System
from pathway import LinearPathway

l1 = symbols('l1', real=True, positive=True)
F, M, g, t = symbols('F, M, g, t')


def make_wire(n):
    """n defines the number of small masses we want to add in the wire"""

    L = symbols('L', positive=True)
    # Creating Reference frames, angles, coordinates, and speeds
    N = ReferenceFrame('N')
    frames = symbols(f'N1:{n + 2}', cls=ReferenceFrame)
    angles = dynamicsymbols(f'alpha1:{n + 2}')
    betas = symbols(f'beta1:{n + 2}')
    qs = dynamicsymbols(f'q1:{2 * (n + 1) + 1}')
    us = dynamicsymbols(f'u1:{2 * (n + 1) + 1}')
    masses = symbols(f'm1:{n + 2}')
    M = symbols('M')

    # Creating the Points
    O = Point('O')
    ps = symbols(f'P1:{n + 2}', cls=Point)
    print(ps)

    # Setting these points in space
    j = 0
    for i in range(n + 1):
        ps[i].set_pos(O, qs[j] * N.x + qs[j + 1] * N.y)
        j += 2

    # Creating Wall and adding it to the system
    wall = RigidBody('W', masscenter=O, frame=N)
    system = System.from_newtonian(wall)

    # Aligning N with system frame
    N.orient_axis(system.frame, 0, system.frame.z)

    # Rotating the frames wrt system frame
    for i, j in zip(frames, angles):
        i.orient_axis(system.frame, j, system.frame.z)

    # Creating the bodies
    bodies = [
        RigidBody(
            f'B{h}', mass=i * j, masscenter=k, frame=l
        ) for h, i, j, k, l in zip(range(1, n + 2), betas, masses, ps, frames)
    ]

    # Adding bodies to the system
    for body in bodies:
        system.add_bodies(body)

    # Assigning velocities to the bodies
    j = 0
    for i in range(n + 1):
        bodies[i].masscenter.set_vel(system.frame, us[j] * N.x + us[j + 1] * N.y)
        bodies[i].frame.set_ang_vel(frames[i], 0)
        j += 2

    # Adding coordinates, speeds, and kdes to the system
    for k in zip(qs, us):
        system.add_coordinates(k[0])
        system.add_speeds(k[1])
        system.add_kdes(k[1] - k[0].diff(t))

    # Applying gravity
    system.apply_uniform_gravity(-g * wall.y)

    # Creating linear pathways between the points and applying force along them
    linear_points = [O] + list(ps)
    print(linear_points)

    linear_paths = []
    for i in range(1, len(linear_points) - 1):
        lin_paths_a = LinearPathway(linear_points[i - 1], linear_points[i])
        lin_paths_b = LinearPathway(linear_points[i + 1], linear_points[i])
        linear_paths.append(lin_paths_a)
        linear_paths.append(lin_paths_b)
    linear_paths.append(LinearPathway(linear_points[-2], linear_points[-1]))
    print(linear_paths)

    for path in linear_paths:
        system.add_loads(path.to_loads(-F)[1])

    # Validating the system
    system.validate_system()

    # Getting the Kane's equations
    system.form_eoms()
    kanel = system.eom_method.to_linearizer()

    op_vals = {
        qs[i]: 0 for i in range(0, 2 * (n + 1), 2)
    } | {
        qs[i]: -(i - i // 2) * l1 for i in range(1, 2 * (n + 1), 2)
    } | {
        F: M * g
    }

    A, B = kanel.linearize(A_and_B=True, op_point=op_vals)
    A = A.subs({l1: L / (n + 1)})

    return A


def make_wire_3d(n, force):
    """n defines the number of small masses we want to add in the wire"""

    L = symbols('L', positive=True)
    # Creating Reference frames, angles, coordinates, and speeds
    N = ReferenceFrame('N')
    frames = symbols(f'N1:{n + 2}', cls=ReferenceFrame)
    angles = dynamicsymbols(f'alpha1:{n + 2}')
    betas = symbols(f'beta1:{n + 1}')
    qs = dynamicsymbols(f'q1:{3 * (n + 1) + 1}')
    us = dynamicsymbols(f'u1:{3 * (n + 1) + 1}')
    masses = symbols(f'm1:{n + 1}')
    M = symbols('M')
    Mf = symbols('Mf')

    # Creating the Points
    O = Point('O')
    ps = symbols(f'P1:{n + 2}', cls=Point)
    print(ps)

    # Setting these points in space
    j = 0
    for i in range(n + 1):
        ps[i].set_pos(O, qs[j] * N.x + qs[j + 1] * N.y + qs[j + 2] * N.z)
        j += 3

    print(qs)

    # Creating Wall and adding it to the system
    wall = RigidBody('W', masscenter=O, frame=N)
    system = System.from_newtonian(wall)

    # Aligning N with system frame
    N.orient_axis(system.frame, 0, system.frame.z)

    # Rotating the frames wrt system frame
    for i, j in zip(frames, angles):
        i.orient_axis(system.frame, j, system.frame.z)

    # Creating the bodies
    bodies = [
        RigidBody(
            f'B{h}', mass=i * j, masscenter=k, frame=l
        ) for h, i, j, k, l in zip(range(1, n + 1), betas, masses, ps[:-1], frames[:-1])
    ]

    last_body = RigidBody(f'B{n + 1}', mass=M, masscenter=ps[-1], frame=frames[-1])
    bodies.append(last_body)
    print(bodies)

    # Adding bodies to the system
    for body in bodies:
        system.add_bodies(body)

    # Assigning velocities to the bodies
    j = 0
    for i in range(n + 1):
        bodies[i].masscenter.set_vel(system.frame, us[j] * N.x + us[j + 1] * N.y + us[j + 2] * N.z)
        bodies[i].frame.set_ang_vel(frames[i], 0)
        j += 3

    # Adding coordinates, speeds, and kdes to the system
    for k in zip(qs, us):
        system.add_coordinates(k[0])
        system.add_speeds(k[1])
        system.add_kdes(k[1] - k[0].diff(t))

    # Applying gravity
    system.apply_uniform_gravity(-g * wall.z)

    # Creating linear pathways between the points and applying force along them
    linear_points = [O] + list(ps)
    print(linear_points)

    linear_paths = []
    for i in range(1, len(linear_points) - 1):
        lin_paths_a = LinearPathway(linear_points[i - 1], linear_points[i])
        lin_paths_b = LinearPathway(linear_points[i + 1], linear_points[i])
        linear_paths.append(lin_paths_a)
        linear_paths.append(lin_paths_b)
    linear_paths.append(LinearPathway(linear_points[-2], linear_points[-1]))
    print(linear_paths)

    for path in linear_paths[:-1]:
        system.add_loads(path.to_loads(-F)[1])

    # Validating the system
    system.validate_system()

    # Getting the Kane's equations
    system.form_eoms()
    kanel = system.eom_method.to_linearizer()

    op_vals = {
        qs[i]: 0 for i in range(0, 3 * (n + 1), 3)
    } | {
        qs[i]: 0 for i in range(1, 3 * (n + 1), 3)
    } | {
        qs[i]: -(i - i // 3) * l1 / 2 for i in range(2, 3 * (n + 1), 3)
    } | {
        F: force
    }

    A, B = kanel.linearize(A_and_B=True, op_point=op_vals)
    A = A.subs({l1: L / (n + 1)})

    return A, system
