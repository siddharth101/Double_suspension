from sympy import zeros, symbols
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
from sympy.physics.mechanics import ReferenceFrame, Point, RigidBody, Particle, System, LinearPathway, Force
from sympy.physics.mechanics import Body, PinJoint, SphericalJoint, PlanarJoint, PrismaticJoint, JointsMethod, inertia
from sympy.physics.mechanics import LinearSpring, LinearDamper
from sympy.physics.mechanics import TorqueActuator
from sympy import atan2
import sympy as smp

def simple_pendulum():

    q1, q2, q3,q4, q5, q6  = dynamicsymbols('q1:7')
    u1, u2, u3, u4, u5, u6 = dynamicsymbols('u1:7')
    alpha, beta, gamma, omega = dynamicsymbols('alpha beta gamma omega')
    
    l, l_0  = symbols('l, l_0', positive=True)
    t = symbols('t')
    m1, g, k1 = symbols('m1, g, k1', positive=True)
    
    N = ReferenceFrame('N')
    O = Point('O')
    
    system = System(frame=N, fixed_point=O)
    P = ReferenceFrame('P')
    P1 = Point('P1')
    P2 = Point('P2')
    
    #P1.set_pos(O, 0*N.x + 0*N.y + 0*N.z)
    
    P2.set_pos(system.fixed_point, (q2+q4)*system.frame.y + (q1+q3)*system.frame.x)



    block = RigidBody('B', mass=m1, masscenter=P2, frame=P)
    
    system.add_bodies(block)

    #N.orient_axis(system.frame, 0, system.frame.z)
    P.orient_axis(system.frame, alpha, system.frame.z)

    # O.set_vel(N, 0)
    # P2.v2pt_theory(P1, N, P)

    block.masscenter.set_vel(system.frame, u1*system.frame.x + u2*system.frame.y)
    block.frame.set_ang_vel(P, 0)
    system.apply_uniform_gravity(-g * system.frame.y)
    #system.add_loads(Force(block, F*block.x))

    linear_path = LinearPathway(P2, system.fixed_point)
    force = k1*(smp.sqrt((q1+q3)**2 + (q2+q4)**2) - l_0)
    linear_path.to_loads(force)
    system.add_loads(linear_path.to_loads(force)[0])
    system.add_loads(linear_path.to_loads(force)[1])

    system.add_coordinates(q1)
    system.add_speeds(u1)

    system.add_kdes(u1 - q1.diff(t))
    system.add_coordinates(q2)
    system.add_speeds(u2)
    system.add_kdes(u2 - q2.diff(t))

    system.validate_system()
    system.form_eoms()
    kanel = system.eom_method.to_linearizer()
    op_vals = {q1:0, q2:-l, l_0: l-m1*g/k1, q3:0, q4:0}
    A, B = kanel.linearize(A_and_B=True, op_point=op_vals)
    
    return A, B, kanel