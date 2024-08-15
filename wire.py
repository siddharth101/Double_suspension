from sympy import zeros, symbols
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
from sympy.physics.mechanics import ReferenceFrame, Point, RigidBody, Particle, System, LinearPathway, Force


l1 = symbols('l1', real=True, positive=True)
F, M, g, t = symbols('F, M, g, t')
def make_wire(n):
    '''n defines the number of small masses we want to add in the wire''' 
    
    # Creating Reference frames, angles, coordinates and speeds
    N = ReferenceFrame('N')
    frames = symbols('N1:{}'.format(n+2), cls=ReferenceFrame)
    angles = dynamicsymbols('alpha1:{}'.format(n+2))
    betas = symbols('beta1:{}'.format(n+2))
    qs = dynamicsymbols('q1:{}'.format(2*(n+1)+1))
    us = dynamicsymbols('u1:{}'.format(2*(n+1)+1))
    
    masses = symbols('m1:{}'.format(n+2))
    M = symbols('M')
   
    #Creating the Points
    O = Point('O')
    ps = symbols('P1:{}'.format(n+2), cls=Point)
    print(ps)

    # Setting these points in space
    j=0
    for i in range(n+1):
        #print(i)
        ps[i].set_pos(O, qs[j]*N.x + qs[j+1]*N.y)
        j+=2
        
    
    # Creating Wall and adding it to the system 
    wall = RigidBody('W', masscenter=O, frame=N)
    system = System.from_newtonian(wall)
   
    # Aligning N with system frame
    N.orient_axis(system.frame, 0, system.frame.z)
    
    # Rotating the frames wrt system frame  
    [i.orient_axis(system.frame, j, system.frame.z) for i,j in list(zip(frames, angles))]

    # Creating the bodies
    bodies = [RigidBody('B{}'.format(h), mass=i*j, masscenter=k, frame=l) for h,i,j,k,l in list(zip(range(1,n+2), betas,masses, ps, frames))]
    
    # Adding bodies to the system
    for i in bodies:
        system.add_bodies(i)
       

    # Assigning velocities to the bodies
    j = 0
    for i in range(n+1):
        #print(frames[i], bodies[i])
        bodies[i].masscenter.set_vel(system.frame, us[j]*N.x + us[j+1]*N.y)
        bodies[i].frame.set_ang_vel(frames[i], 0)
        j+=2
        
    # Adding co-ordinates, speeds and kdes to the system 
    for k in list(zip(qs, us)):
        system.add_coordinates(k[0])
        system.add_speeds(k[1])
        system.add_kdes(k[1]-k[0].diff(t))
       
    # Applying gravity
    system.apply_uniform_gravity(-g * wall.y)
   
    # Now we will create linearpathways between the points and apply force along them
    linear_points = [O] + [i for i in ps]
    print(linear_points)
    
    linear_paths = []
    for i in range(1, len(linear_points)-1):
        lin_paths_a = LinearPathway(linear_points[i-1], linear_points[i])
        lin_paths_b = LinearPathway(linear_points[i+1], linear_points[i])
        linear_paths.append(lin_paths_a)
        linear_paths.append(lin_paths_b)
    linear_paths.append(LinearPathway(linear_points[-2], linear_points[-1]))
    print(linear_paths)
    
    [system.add_loads(i.to_loads(-F)[1]) for i in linear_paths]
   

    # System is prepared, we can validate it
    system.validate_system()
    
    # Getting the kane's equations
    system.form_eoms()
    kanel = system.eom_method.to_linearizer()
    
    op_vals = {qs[i]:0 for i in range(0,2*(n+1),2)} |  {qs[i]:-(i-i//2)*l1 for i in range(1,2*(n+1),2)} | {F:M*g}
    
    A, B = kanel.linearize(A_and_B=True, op_point=op_vals)
    
    return A
