from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
from sympy.physics.mechanics import ReferenceFrame, Point, RigidBody
from sympy import zeros, symbols

import numpy as np
import sympy as sp
import sympy.physics.mechanics as me
from IPython.display import display
from suspycious import Model
import suspycious.components as scmp

k1, k2, k3, k4,  g, c1, c2, t = symbols('k1 k2 k3 k4 g c1 c2 t')

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
D = model2.add(scmp.RigidBody("D"))
F = model2.add(scmp.RigidBody("F"))


l1, l2, l3, l4 = symbols('l1 l2 l3 l4')

S.set_global_position(0,0,0)
B.set_global_position(0,0, -l1)
C.set_global_position(0,0, -l1 -l2)
D.set_global_position(0,0, -l1 -l2 -l3)
F.set_global_position(0,0, -l1 -l2 -l3 -l4)

S.set_global_orientation(0,0,0)
B.set_global_orientation(0,0,0)
C.set_global_orientation(0,0,0)
D.set_global_orientation(0,0,0)
F.set_global_orientation(0,0,0)




#Define geometrical properties (see images)
d1, d2, n1, n2, w1, w2 = symbols('d1 d2 n1 n2 w1 w2', real=True, positive=True)
d3, n3, w3 = symbols('d3 n3 w3', real=True, positive=True)
d4, n4, w4 = symbols('d4 n4 w4', real=True, positive=True)
d5, n5, w5 = symbols('d5 n5 w5', real=True, positive=True)
d6, n6, w6 = symbols('d6 n6 w6', real=True, positive=True)
d7, n7, w7 = symbols('d7 n7 w7', real=True, positive=True)
d8, d9 = symbols('d8 d9', real=True, positive=True)


def add_points(body, point='P', attachment_points=[w1, n1, d1]):
    points = [Point('{}{}'.format(point,j)) for j in range(1,5)]
    
    points_ = [str(i) for i in points]
    
    w_ = attachment_points[0]
    n_ = attachment_points[1]
    d_ = attachment_points[2]
    
    
    body.add_fixed_point(points_[0], dx = w_, dy = -n_, dz = d_)
    body.add_fixed_point(points_[1], dx = w_, dy = n_, dz = d_)
    body.add_fixed_point(points_[2], dx = -w_, dy = -n_, dz = d_)
    body.add_fixed_point(points_[3], dx = -w_, dy = n_, dz = d_)

    return

def get_tension_dirs(body1, body2, points_1='top', points_2='top'):
    
    body1_points = list(body1.points.values())
    body2_points = list(body2.points.values())
    
    if points_1 == 'top':
        body1_points = body1_points[:4]
    else:
        body1_points = body1_points[4:]
    if points_2 == 'top':
        body2_points = body2_points[:4]
    else:
        body2_points = body2_points[4:]
        
    print(body1_points)
    print(body2_points)
    
    ten_dirs = [i.pos_from(j).express(S.global_frame).normalize() for i,j in list(zip(body1_points, body2_points))]
        
    return ten_dirs
    
    
def give_deltas(body1, body2, points_1='top', points_2='top'):
    
    body1_points = list(body1.points.values())
    body2_points = list(body2.points.values())
    
    if points_1 == 'top':
        body1_points = body1_points[:4]
    else:
        body1_points = body1_points[4:]
    if points_2 == 'top':
        body2_points = body2_points[:4]
    else:
        body2_points = body2_points[4:]
        
    points_ = list(zip(body1_points, body2_points))
    
    print(points_)
    
    deltas_ = [i.pos_from(j).express(S.global_frame).magnitude() - i.pos_from(j).express(S.global_frame).subs(model2.op_point).magnitude() 
              for i,j in points_]
    return deltas_
    
    
def give_tensions(n_body=1, k=k1, delta_values=None):
    
    n = n_body -1
    
    #masses = [B.M, C.M, D.M, F.M]
    masses_ = masses[n:]
    
    tensions = [sp.Rational(1,4)*(np.sum(masses_))*g + k*i for i in delta_values]
        
    return tensions 



add_points(body=S, point='P', attachment_points=[w1, n1, 0])

add_points(body=B, point='B', attachment_points=[w1, n1, d1])
add_points(body=B, point='C', attachment_points=[w2, n2, -d2])

add_points(body=C, point='D', attachment_points=[w2, n2, d3])
add_points(body=C, point='E', attachment_points=[w3, n3, -d4])

add_points(body=D, point='F', attachment_points=[w3, n3, d5])
add_points(body=D, point='G', attachment_points=[w4, n4, -d6])

add_points(body=F, point='H', attachment_points=[w4, n4, d7])
add_points(body=F, point='J', attachment_points=[w5, n5, -d8])


dirT11, dirT12, dirT13, dirT14 = get_tension_dirs(body1=S,points_1='top', body2=B, points_2='top')
dirT21, dirT22, dirT23, dirT24 = get_tension_dirs(body1=B,points_1='bottom', body2=C, points_2='top')

dirT31, dirT32, dirT33, dirT34 = get_tension_dirs(body1=C,points_1='bottom', body2=D, points_2='top')
dirT41, dirT42, dirT43, dirT44 = get_tension_dirs(body1=D,points_1='bottom', body2=F, points_2='top')


delta_SB = give_deltas(body1=S, body2=B )
delta_BC = give_deltas(body1=B, body2=C, points_1='bottom', points_2='top' )
delta_CD = give_deltas(body1=C, body2=D, points_1='bottom', points_2='top' )
delta_DF = give_deltas(body1=D, body2=F, points_1='bottom', points_2='top' )

masses = [B.M, C.M, D.M, F.M]


T11, T12, T13, T14 = give_tensions(n_body=1, k=k1, delta_values=delta_SB)
T21, T22, T23, T24 = give_tensions(n_body=2, k=k2, delta_values=delta_BC)
T31, T32, T33, T34 = give_tensions(n_body=3, k=k3, delta_values=delta_CD)
T41, T42, T43, T44 = give_tensions(n_body=4, k=k4, delta_values=delta_DF)

print("Prepared the system, now applying forces")

B.body.apply_force(T11*dirT11 , point=B.B1)
B.body.apply_force(T12*dirT12 , point=B.B2)
B.body.apply_force(T13*dirT13 , point=B.B3)
B.body.apply_force(T14*dirT14 , point=B.B4)
B.body.apply_force(-T21*dirT21 , point=B.C1)
B.body.apply_force(-T22*dirT22 , point=B.C2)
B.body.apply_force(-T23*dirT23 , point=B.C3)
B.body.apply_force(-T24*dirT24 , point=B.C4)
B.Fz = -np.sum([i*g for i in masses])

C.body.apply_force(T21*dirT21 , point=C.D1)
C.body.apply_force(T22*dirT22 , point=C.D2)
C.body.apply_force(T23*dirT23 , point=C.D3)
C.body.apply_force(T24*dirT24 , point=C.D4)
C.body.apply_force(-T31*dirT31 , point=C.E1)
C.body.apply_force(-T32*dirT32 , point=C.E2)
C.body.apply_force(-T33*dirT33 , point=C.E3)
C.body.apply_force(-T34*dirT34 , point=C.E4)
C.Fz = -np.sum([i*g for i in masses[1:]])

D.body.apply_force(T31*dirT31 , point=D.F1)
D.body.apply_force(T32*dirT32 , point=D.F2)
D.body.apply_force(T33*dirT33 , point=D.F3)
D.body.apply_force(T34*dirT34 , point=D.F4)
D.body.apply_force(-T41*dirT41 , point=D.G1)
D.body.apply_force(-T42*dirT42 , point=D.G2)
D.body.apply_force(-T43*dirT43 , point=D.G3)
D.body.apply_force(-T44*dirT44 , point=D.G4)
D.Fz = -np.sum([i*g for i in masses[2:]])

F.body.apply_force(T41*dirT41 , point=F.H1)
F.body.apply_force(T42*dirT42 , point=F.H2)
F.body.apply_force(T43*dirT43 , point=F.H3)
F.body.apply_force(T44*dirT44 , point=F.H4)
F.Fz = -np.sum([i*g for i in masses[3:]])


print("Extracting State Space")
import time
tic = time.time()
kane = model2.extract_statespace()
tac = time.time()
print(tac - tic)
print("Done")


A, B = kane.A, kane.B

L1, L2, L3, L4 = symbols('L1, L2, L3, L4', positive=True)

A = A.subs({-d1 + l1:L1, -d2-d3+l2:L2, -d4 -d5 +l3:L3, -d6 -d7 + l4:L4})
B = B.subs({-d1 + l1:L1, -d2-d3+l2:L2, -d4 -d5 +l3:L3, -d6 -d7 + l4:L4})

print(A[40,10])