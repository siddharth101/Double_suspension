from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
from sympy import zeros, symbols

import numpy as np
import sympy as sp
import sympy.physics.mechanics as me
from IPython.display import display
from suspycious import Model
import suspycious.components as scmp

from sympy import atan2
import sympy as smp
from sympy import zeros, symbols
from sympy import symbols, Matrix, solve, simplify
from sympy import Matrix
from sympy.physics.mechanics import Body, PinJoint, SphericalJoint, PlanarJoint, PrismaticJoint, JointsMethod, inertia
from sympy.physics.mechanics import dynamicsymbols
from sympy import Symbol
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, outer
from sympy.physics.mechanics import ReferenceFrame, Point, RigidBody
from sympy.physics.mechanics import kinetic_energy, potential_energy, Point, Particle


from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
from sympy import zeros, symbols

import numpy as np
import sympy as sp
import sympy.physics.mechanics as me
from IPython.display import display
from suspycious import Model
import suspycious.components as scmp

model2 = Model()
S = model2.add(scmp.RigidBody("S"))
B = model2.add(scmp.RigidBody("B"))
C = model2.add(scmp.RigidBody("C"))

q1, q2, q3 = dynamicsymbols('q1 q2 q3', positive=True)
w1, w2, d1, d2, n1, n2, l1, l2 = symbols('w1 w2 d1 d2 n1 n2 l1 l2', real=True, positive=True)

S.set_global_position(0,0,0)
B.set_global_position(0,0, 0)
#C.set_global_position(0,0, -l1)
C.set_global_position(0, 0, -l1)

S.set_global_orientation(0,0,0)
B.set_global_orientation(0,0,0)
C.set_global_orientation(0,0,0)

alpha_z1, beta_z1 = dynamicsymbols('alpha_z1 beta_z1')
alpha_x1, beta_x1 = dynamicsymbols('alpha_x1 beta_x1')
alpha_y1, beta_y1 = dynamicsymbols('alpha_y1 beta_y1')


omega_z1 = dynamicsymbols('omega_z1')
omega_x1 = dynamicsymbols('omega_x1')
omega_y1 = dynamicsymbols('omega_y1')
q2, q3, q4 = dynamicsymbols('q2, q3, q4', positive=True)
u2, u3, u4 = dynamicsymbols('u2 u3, u4')

T1, k1, delta_l1, x_0, t = symbols('T1, k1, delta_l1, x_0, t')


Az = C.com.pos_from(S.com).dot(S.global_frame.z)
Ax = C.com.pos_from(S.com).dot(S.global_frame.x)
Ay = C.com.pos_from(S.com).dot(S.global_frame.y)

r = C.com.pos_from(S.com).magnitude()

B.frame.dcm(S.frame)


#angle_springx = atan2(smp.sqrt(r**2 - Az**2), Az)
angle_springx = atan2(Ay, Az)
angle_springz = atan2(Ax, Ay)
rev1 = SphericalJoint(name='p1p2', parent=S.body, child=B.body, parent_point=S.com, child_point=B.com, amounts=[angle_springx,0,0],
                      coordinates=[alpha_x1,alpha_y1, alpha_z1],
                      speeds=[omega_x1, omega_y1, omega_z1],
                     rot_order='XYZ')

B.frame.dcm(S.frame)
C.frame.dcm(S.frame)

C.frame.dcm(B.frame)
pos_c_b = C.com.pos_from(B.com).express(model2.frame).subs(S.op_point).subs(B.op_point).magnitude().subs({C.dy:0, C.dz:0, C.dx:0, C.x:0})

C.frame.orient_axis(B.frame, 0, B.frame.z)
C.frame.orient_body_fixed(B.frame, angles=(0,0,0), rotation_order='XYZ')
rev2 = PrismaticJoint('J2', parent=B.body, child=C.body, parent_point=B.com, child_point=C.com, coordinates=q4, speeds=u4,
                    joint_axis=-C.frame.y)

g = symbols('g')
C.body.apply_force(-C.body.mass*g*S.global_frame.z, point=C.com, reaction_body=B.body, reaction_point=B.com)
T1 = k1*delta_l1.subs({delta_l1:q4 - x_0})#.subs({q4:{smp.sqrt(C.y**2 + C.z**2)}})

T1 = T1.subs({q4:pos_c_b})

C.body.apply_force(T1*C.frame.z, point=C.com, reaction_body=B.body, reaction_point=B.com)
method = JointsMethod(B.body, rev2)
method.form_eoms()


kane = model2.extract_statespace()
kaneeq = kane.kane.kanes_equations()


A = kane.A
B = kane.B

