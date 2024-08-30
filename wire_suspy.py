from sympy import zeros, symbols
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
from sympy.physics.mechanics import ReferenceFrame, Point, RigidBody, Particle, System, Force
from pathway import LinearPathway
import numpy as np


def make_wire_suspy(n, linpath, index, global_frame, suspension_body):
    
    point1 = linpath['paths'][index].attachments[0]
    point2 = linpath['paths'][index].attachments[1]
    
    print(point1, point2)
    
    bodies_str = ''.join([point1.name[0], point2.name[0]])
    
    body1 = linpath['bodies'][0]
    body2 = linpath['bodies'][-1]
    
    mass1 = linpath['bodies'][0].mass
    mass1 = linpath['bodies'][1].mass
   
    gb_frame = global_frame
    ind_str = index+1
    l1 = symbols('l1', positive=True)
    L = symbols('L_' + bodies_str+ '{}{}'.format(index+1, index+1), positive=True)
    F,g,t = symbols('F, g, t')
    
    
    frames = symbols('N{}{}_{}1:{}'.format(ind_str, ind_str,bodies_str,n+2), cls=ReferenceFrame)
    angles = dynamicsymbols('alpha{}{}_{}1:{}'.format(ind_str, ind_str,bodies_str,n+2))
    betas = symbols('beta{}{}_{}1:{}'.format(ind_str, ind_str,bodies_str,n+2))
    #qs = dynamicsymbols('q1:{}'.format(3*(n+1)+1))
    qs = dynamicsymbols('q{}{}_{}1:{}'.format(ind_str, ind_str,bodies_str, 3*(n+1)+1))
    us = dynamicsymbols('u{}{}_{}1:{}'.format(ind_str, ind_str,bodies_str, 3*(n+1)+1))
    print(qs)
    masses = symbols('m{}{}_{}1:{}'.format(ind_str, ind_str, bodies_str, n+2))
    M = symbols('M')
    
    
    model = suspension_body.model
    
    op_vals = {qs[i]:0 for i in range(0,3*(n+1),3)} | \
    {qs[i]:0 for i in range(1,3*(n+1),3)} | \
    {qs[i]:-(i-i//3)*l1/2 for i in range(2,3*(n+1),3)} | \
    {F:mass1*g} 
    
    
    dist21 = point2.pos_from(point1).express(gb_frame)
    distx = dist21.dot(gb_frame.x)
    disty = dist21.dot(gb_frame.y)
    distz = dist21.dot(gb_frame.z)
    
    angle = distz/dist21.magnitude()
    angles[-1] = angle - np.sum([i for i in angles[:-1]])
    
    
    S = suspension_body
    wall = RigidBody('W', masscenter=S.global_origin, frame=gb_frame)
    system = System.from_newtonian(wall)
    
    #print(angles)
    [i.orient_axis(system.frame, j, system.frame.z) for i,j in list(zip(frames, angles))]
    
    
    ps = symbols('P{}{}_{}1:{}'.format(ind_str, ind_str, bodies_str, n+2), cls=Point)



    # Setting these points in space
    j=0
    for i in range(n+1):
        print(ps[i], qs[j], qs[j+1], qs[j+2])
        ps[i].set_pos(point1, qs[j]*gb_frame.x + qs[j+1]*gb_frame.y + qs[j+2]*gb_frame.z)
        j+=3

    
    op_point_wire = {qs[-3]:distx, qs[-2]:disty, qs[-1]:distz}
    
    
    # Creating the bodies
    bodies = [RigidBody('B{}{}_{}{}'.format(ind_str, ind_str,bodies_str,h), mass=i*j, masscenter=k, frame=l)
              for h,i,j,k,l in list(zip(range(1,n+2), betas,masses, ps, frames))]
    
 
    # Adding bodies to the system
    for i in bodies:
        system.add_bodies(i)
        
    # Assigning velocities to the bodies
    j = 0
    for i in range(n+1):
        print(frames[i], bodies[i])
        bodies[i].masscenter.set_vel(system.frame, 
                                     us[j]*gb_frame.x + us[j+1]*gb_frame.y + us[j+2]*gb_frame.z)
        bodies[i].frame.set_ang_vel(frames[i], 0)
        j+=3
        
        

    #print(list(zip(qs, us)))
    for k in list(zip(qs, us)):
        #print(qs)
        system.add_coordinates(k[0])
        system.add_speeds(k[1])
        system.add_kdes(k[1]-k[0].diff(t))
        
    # Applying gravity
    system.apply_uniform_gravity(-g * wall.z)
    
   
    
    print(ps[-1].pos_from(point2).express(gb_frame).subs({qs[-3]:distx, qs[-2]:disty, qs[-1]:distz}).subs(model.op_point))
    

    
    linear_points = [point1] + [i for i in ps]# + [point2]
    print(linear_points)
    
    linear_paths = []
    for i in range(1, len(linear_points)-1):
        lin_paths_a = LinearPathway(linear_points[i-1], linear_points[i])
        lin_paths_b = LinearPathway(linear_points[i+1], linear_points[i])
        linear_paths.append(lin_paths_a)
        linear_paths.append(lin_paths_b)
    linear_paths.append(LinearPathway(linear_points[-2], linear_points[-1]))
    print(linear_paths)
    [system.add_loads(i.to_loads(-F)[1]) for i in linear_paths[:-1]]
    
    system.validate_system()
    
    # Getting the kane's equations
    system.form_eoms()
    kanel = system.eom_method.to_linearizer()
    
    
    A, B = kanel.linearize(A_and_B=True, op_point=op_vals)
    
    A = A.subs({l1: L/(n+1)}) 
    B = B.subs({l1:L/(n+1)}) 
    
    
    return system, A
