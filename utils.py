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

def get_tension_dirs(body1, body2, suspension_body, points_1='top', points_2='top'):
    
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
    
    ten_dirs = [i.pos_from(j).express(suspension_body.global_frame).normalize() for i,j in list(zip(body1_points, body2_points))]
        
    return ten_dirs
    
    
def give_deltas(body1, body2, suspension_body,model,points_1='top', points_2='top'):
    
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
    
    op_point = model.op_point
    sb = suspension_body
    print(points_)
    
    deltas_ = [i.pos_from(j).express(sb.global_frame).magnitude() - i.pos_from(j).express(sb.global_frame).subs(op_point).magnitude() 
              for i,j in points_]
    return deltas_
    
    
def give_tensions(masses, n_body=1, k=k1, delta_values=None ):
    
    n = n_body -1
    
    #masses = [B.M, C.M, D.M, F.M]
    masses_ = masses[n:]
    
    tensions = [sp.Rational(1,4)*(np.sum(masses_))*g + k*i for i in delta_values]
        
    return tensions 
