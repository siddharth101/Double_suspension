from sympy import zeros, symbols
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
from sympy.physics.mechanics import ReferenceFrame, Point, RigidBody, Particle, System, LinearPathway, Force
from sympy.physics.mechanics import Body, PinJoint, SphericalJoint, PlanarJoint, PrismaticJoint, JointsMethod, inertia
from sympy.physics.mechanics import LinearSpring, LinearDamper
from sympy.physics.mechanics import TorqueActuator
from sympy import atan2
import sympy as smp
import numpy as np

def particle_pend(n, force=None):
    """n defines the number of small masses we want to add in the wire"""

    L = symbols('L', positive=True)
    # Creating Reference frames, angles, coordinates, and speeds
    N = ReferenceFrame('N')
    frames = symbols(f'N1:{n + 2}', cls=ReferenceFrame)
    angles = dynamicsymbols(f'alpha1:{n + 2}')
    betas = symbols(f'beta1:{n + 2}')
    qs = dynamicsymbols(f'q1:{2 * (n + 1) + 1}')
    us = dynamicsymbols(f'u1:{2 * (n + 1) + 1}')
    qc = dynamicsymbols(f'qc1:4')
    uc = dynamicsymbols(f'uc1:4')
    masses = symbols(f'm1:{n + 1}')
    M = symbols('M')
    masses = list(masses) + [M]
    
    metric = [1,-1,0]
   
    t,g = symbols('t,g')
    l1, l2 = symbols('l1 l2', positive=True)
    l_s = symbols('l_1:{}'.format(n+2))
    ks = symbols('k1:{}'.format(n+2))
    
    # Creating the Points
    O = Point('O')
    ps = symbols(f'P1:{n + 2}', cls=Point)
   # print(ps)

    system = System(fixed_point=O, frame=N)
    
    # Setting these points in space
    j = 0
    for i in range(n):
        ps[i].set_pos(system.fixed_point, (qs[j]) * metric[0] * system.frame.x +
                      (qs[j + 1]) * metric[1] * system.frame.y)
        j += 2

    ps[-1].set_pos(system.fixed_point, (qs[-2] + qc[0]) * metric[0] * system.frame.x +
                      (qs[-1] + qc[1]) * metric[1] * system.frame.y)
    
    
    # Creating Wall and adding it to the system
    #wall = RigidBody('W', masscenter=O, frame=N)
    #system = System.from_newtonian(wall)

    
    
    # Aligning N with system frame
    #N.orient_axis(system.frame, 0, system.frame.z)

    # Rotating the frames wrt system frame
    for i, j in zip(frames, angles):
        i.orient_axis(system.frame, j, system.frame.z)
    
    bodies = [
        Particle(
            f'B{h}', mass=i, point=j
        ) for h, i, j in zip(range(1, n + 1), masses, ps[:-1])
    ]
    
    
    last_body = RigidBody(
            f'B{n+1}', mass=masses[-1], masscenter=ps[-1], frame=frames[-1])
    
    bodies = bodies + [last_body]
    # Adding bodies to the system
    for body in bodies:
        system.add_bodies(body)
    #system.add_bodies(last_body)

    #Assigning velocities to the bodies
    j = 0
    for i in range(n + 1):
        bodies[i].masscenter.set_vel(system.frame, (us[j]) * system.frame.x + (us[j + 1]) * system.frame.y)
        #bodies[i].frame.set_ang_vel(frames[i], 0)
        j += 2
    

    
    bodies[-1].masscenter.set_vel(system.frame, (us[-2]) * system.frame.x + (us[-1]) * system.frame.y)
    bodies[-1].frame.set_ang_vel(frames[-1], 0)


    # Adding coordinates, speeds, and kdes to the system
    for k in zip(qs, us):
        system.add_coordinates(k[0])
        system.add_speeds(k[1])
        system.add_kdes(k[1] - k[0].diff(t))

    # Applying gravity
    system.apply_uniform_gravity(-g * system.frame.y)

    # Creating linear pathways between the points and applying force along them
    linear_points = [O] + list(ps)
    #print(linear_points)

    linear_paths = []
    for i in range(1, len(linear_points) - 1):
        lin_paths_a = LinearPathway(linear_points[i - 1], linear_points[i])
        lin_paths_b = LinearPathway(linear_points[i + 1], linear_points[i])
        linear_paths.append(lin_paths_a)
        linear_paths.append(lin_paths_b)
    linear_paths.append(LinearPathway(linear_points[-2], linear_points[-1]))
    #print(linear_paths)

    
    deltas_ = []
    delta_iterator = 0
    for paths in linear_paths[::2]:
        delta_val = paths.attachments[0].pos_from(paths.attachments[1]).magnitude() - l_s[delta_iterator]
        deltas_.append(delta_val)
        delta_iterator+=1
    
    forces_values = [i*j for i,j in list(zip(list(ks), deltas_))]
    
    
    F = symbols('F')
    #forces_ = symbols('{}1:{}'.format(force,n+2))
    
    
    forces_ = forces_values
    forces_a = [forces_[0]] + [val for val in forces_[1:] for _ in (0,1)]
    #print(forces_)
    force_iterator = 0
    for path in linear_paths:
        system.add_loads(path.to_loads(-forces_a[force_iterator])[1])
        force_iterator+=1
        #system.add_loads(path.to_loads(-F)[1])

    # Validating the system
    system.validate_system()

    # Getting the Kane's equations
    system.form_eoms()
    kanel = system.eom_method.to_linearizer()

    op_vals = {
        qs[i]: 0 for i in range(0, 2 * (n + 1), 2)
    } | {
        qs[i]: (i - i // 2) * l1 for i in range(1, 2 * (n + 1), 2)
    } | {
        qc[i]:0 for i in range(0,3,1)
    } 
    
    # |  {
    #     F: force
    # }
    

    A, B = kanel.linearize(A_and_B=True, op_point=op_vals)
    
    # print([np.sum(masses[i:]) for i in range(n+1)])
  
    lsdict = {i:l1 - (np.sum(masses[m:])*g)/k for i,m,k in list(zip(list(l_s),  range(n+2), list(ks))) }
   # print(lsdict)
    
    A = A.subs(lsdict)
    B = B.subs(lsdict)
    
    #print(kanel.r)
    
    A = A.subs({l1: L / (n + 1)})

    return A,B, system
