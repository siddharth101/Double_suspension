from sympy.physics.mechanics import Point, RigidBody
from sympy import zeros, symbols, Matrix
import suspycious.components as scmp
from pathway import LinearPathway
import sympy as sp
import numpy as np


def get_tension_dirs(body1, body2, points_1='top', points_2='top', suspension_body=None):
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
    S = suspension_body
    ten_dirs = [
        i.pos_from(j).express(S.global_frame).normalize()
        for i, j in zip(body1_points, body2_points)
    ]

    return ten_dirs


def get_wire_dir(body1, body2, points_1='top', points_2='top', suspension_body=None):
    paths = get_linpaths(
        body1, body2, points_1=points_1, points_2=points_2,
        suspension_body=suspension_body
    )['paths']

    S = suspension_body

    dirs = [
        i.attachments[0].pos_from(i.attachments[1]).express(S.global_frame).normalize()
        for i in paths
    ]

    return dirs


# def get_linpaths(body1, body2, points_1='top', points_2='top', suspension_body=None):
#     body1_points = list(body1.points.values())
#     body2_points = list(body2.points.values())

#     if points_1 == 'top':
#         body1_points = body1_points[:4]
#     else:
#         body1_points = body1_points[4:]
#     if points_2 == 'top':
#         body2_points = body2_points[:4]
#     else:
#         body2_points = body2_points[4:]

#     print(body1_points)
#     print(body2_points)
#     S = suspension_body
#     linpaths = {
#         'paths': [LinearPathway(i, j) for i, j in zip(body1_points, body2_points)],
#         'bodies': [body1.body, body2.body]
#     }

#     return linpaths


def give_deltas(body1, body2, suspension_body, model, points_1='top', points_2='top'):
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
    deltas_ = [
        i.pos_from(j).express(sb.global_frame).magnitude() -
        i.pos_from(j).express(sb.global_frame).subs(op_point).magnitude()
        for i, j in points_
    ]
    return deltas_


def add_attachment_points(name, body, position='top', *args, **kwargs):
    points = [Point(f'{name}{j}') for j in range(1, 9)]
    points_ = [str(i) for i in points]  # these are the points
    matrix_attachment = Matrix([(1, -1, 1), (1, 1, 1), (-1, -1, 1), (-1, 1, 1),
                                (1, -1, -1), (1, 1, -1), (-1, -1, -1), (-1, 1, -1)])

    if position == 'top':
        for i in range(4):
            body._unprotect()
            pt = points_[i]
            dx_pt = body._add_symbol(f"{pt}_x", f"{pt}_x", level=0, real=True)
            dy_pt = body._add_symbol(f"{pt}_y", f"{pt}_y", level=0, real=True)
            dz_pt = body._add_symbol(f"{pt}_z", f"{pt}_z", level=0, real=True)
            attachment_points = [dx_pt, dy_pt, dz_pt]
            matrix_col_1 = matrix_attachment[i:i+1, 0] * attachment_points[0]
            matrix_col_2 = matrix_attachment[i:i+1, 1] * attachment_points[1]
            matrix_col_3 = matrix_attachment[i:i+1, 2] * attachment_points[2]
            a = [matrix_col_1[0, 0], matrix_col_2[0, 0], matrix_col_3[0, 0]]
            print(a[0], a[1], a[2])
            body.add_fixed_point(points_[i], a[0], a[1], a[2])

    else:
        for i in range(4, 8):
            body._unprotect()
            pt = points_[i]
            dx_pt = body._add_symbol(f"{pt}_x", f"{pt}_x", level=0, real=True)
            dy_pt = body._add_symbol(f"{pt}_y", f"{pt}_y", level=0, real=True)
            dz_pt = body._add_symbol(f"{pt}_z", f"{pt}_z", level=0, real=True)
            attachment_points = [dx_pt, dy_pt, dz_pt]

            matrix_col_1 = matrix_attachment[i:i+1, 0] * attachment_points[0]
            matrix_col_2 = matrix_attachment[i:i+1, 1] * attachment_points[1]
            matrix_col_3 = matrix_attachment[i:i+1, 2] * attachment_points[2]
            a = [matrix_col_1[0, 0], matrix_col_2[0, 0], matrix_col_3[0, 0]]
            print(a[0], a[1], a[2])
            body.add_fixed_point(points_[i], a[0], a[1], a[2])

    body._protect()


def give_tensions(n_body, k, masses, delta_values=None):
    n = n_body - 1
    g = symbols('g')
    masses_ = masses[n:]

    tensions = [sp.Rational(1, 4) * (np.sum(masses_)) * g + k * i for i in delta_values]

    return tensions


def add_points(body, point, attachment_points):
    points = [Point(f'{point}{j}') for j in range(1, 5)]
    points_ = [str(i) for i in points]

    w_ = attachment_points[0]
    n_ = attachment_points[1]
    d_ = attachment_points[2]

    body.add_fixed_point(points_[0], dx=w_, dy=-n_, dz=d_)
    body.add_fixed_point(points_[1], dx=w_, dy=n_, dz=d_)
    body.add_fixed_point(points_[2], dx=-w_, dy=-n_, dz=d_)
    body.add_fixed_point(points_[3], dx=-w_, dy=n_, dz=d_)

    return


def apply_force(body, forcedict, global_frame, bottom=True):
    top_path = forcedict[body]['paths'][0]['paths']
    top_points = list(forcedict[body]['points'].values())[:4]
    top_tensions = forcedict[body]['tensions'][:4]

    for i in range(4):
        body.body.apply_force(
            top_path[i].to_loads(-top_tensions[i], frame=global_frame)[1][1],
            point=top_points[i]
        )

    if bottom:
        bottom_path = forcedict[body]['paths'][1]['paths']
        bottom_points = list(forcedict[body]['points'].values())[4:]
        bottom_tension = forcedict[body]['tensions'][4:]

        for i in range(4):
            body.body.apply_force(
                bottom_path[i].to_loads(-bottom_tension[i], frame=global_frame)[0][1],
                point=bottom_points[i]
            )
    else:
        pass

    return


def make_wire_bodies(n, model, top_body, bottom_body, attachment_points, suspension_body):

    
    global_frame = suspension_body.global_frame
    global_origin = suspension_body.global_origin
    
    name_ = top_body.body.name[0] + bottom_body.body.name[0]
    num = 4*n
    
    names = [name_+'{}{}'.format(j+1,i) for j in range(n) for i in range(1,4+1)]
    print(names)
    for name in names:
        model.add(scmp.RigidBody(name))
    
    print("Added {} bodies to the model between {} and {}".format(num, top_body.body.name, bottom_body.body.name))
    
    top_frame = top_body.frame
    bottom_frame = bottom_body.frame
    
    bottom_body_bottom_points = list(bottom_body.points)[:4]
    
    p1, p2 = attachment_points
    
    d1_ = p1.pos_from(global_origin).dot(global_frame.z).subs(model.op_point)
    print(d1_)
    d2_ = p2.pos_from(p1).dot(global_frame.z).subs(model.op_point)
    print(d2_)
    points_wire_coeff = [sp.Rational(1+i,n+1) for i in range(n)]
    
    d_ = [d1_ + i*d2_ for i in points_wire_coeff]
    print(d_)
    
    coeffs = []
    for point in bottom_body_bottom_points:
        coeff_ = [bottom_body.points[point].pos_from(bottom_body.body.masscenter).subs(model.op_point).dot(bottom_frame.x),
        bottom_body.points[point].pos_from(bottom_body.body.masscenter).subs(model.op_point).dot(bottom_frame.y)]
    
        coeffs.append(coeff_)
        
    coeffs_new = []
    for i in range(len(d_)):
        coeffs_new.append([j + [d_[i]]  for j in coeffs])
    
    
    coeffs_new = [item for items in coeffs_new for item in items]
    
    
    added_bodies = [i[-1] for i in list(model.components.items())[-num:]]
    added_coms = [i[-1].com for i in list(model.components.items())[-num:]]
    #print(added_bodies)
    
    ### Setting global position and global orientation 
    for i in range(num):
        added_bodies[i].set_global_position(coeffs_new[i][0], coeffs_new[i][1], coeffs_new[i][2])
        added_bodies[i].set_global_orientation(0,0,0)
    
    if top_body!=suspension_body:
        top_points = list(top_body.points)[4:]
        top_points_ = [top_body.points[i] for i in top_points]
    else:
        top_points = list(top_body.points)
        top_points_ = [top_body.points[i] for i in top_points]
    
    bottom_points = list(bottom_body.points)[:4]
    bottom_points_ = [bottom_body.points[i] for i in bottom_points]
    
    print(added_coms)
    points_list = top_points_ + added_coms + bottom_points_
    print(added_bodies)
    
    body_list = 4*[top_body] + added_bodies + 4*[bottom_body]
    #return coeffs_new    
    
    points_body_dict = {'points':points_list, 'bodies':body_list}
    
    return points_body_dict
    
def get_points(n, model, num_blocks, bodies, attach_points, suspension_body):

    points_bodies = {}

    #bodies = [[S, B], [B,C], [C,D], [D,F]]
    #attach_points = [[S.P1, B.B1], [B.C1, C.D1], [C.E1, D.F1], [D.G1, F.H1]]
    for i in range(num_blocks):
        points_bodies['body{}'.format(i+1)] = make_wire_bodies(n=n, model=model, top_body=bodies[i][0], bottom_body=bodies[i][1],
                                                               attachment_points=[attach_points[i][0], attach_points[i][1]], suspension_body=suspension_body)

    return points_bodies

def get_deltas(points_body, model, suspension_body):
    
    n = len(points_body['points'])//4 - 1
    
    deltas_ = []
    op_point = model.op_point
    sb = suspension_body
    deltas_a = {}
    top_body_name = points_body['bodies'][0].name
    bottom_body_name = points_body['bodies'][-1].name
    delta_body_str = [f'del_{top_body_name}_{bottom_body_name}_{i+1}' for i in range(len(points_body['bodies'])-4)]

    for j in range(4): #4 because of 4 wires
        for i in range(n): #n is wire parts per wire (Depends on violin modes)
            m = 4*i + j
            attach_points = (points_body['points'][m], points_body['points'][m+4])
            print(attach_points)
            delta = [
        attach_points[0].pos_from(attach_points[1]).express(sb.global_frame).magnitude() -
        attach_points[0].pos_from(attach_points[1]).express(sb.global_frame).subs(op_point).magnitude()]    

            deltas_.append(delta)
    

    deltas_unwrapped = [item for items in deltas_ for item in items]

    dict_del_body = dict(list(zip(delta_body_str, deltas_unwrapped)))
        
    return  dict_del_body#deltas_unwrapped#dict_del_body #deltas_unwrapped
    
    
def deltas_dict(points_body, num_blocks, model, suspension_body):

    delta_ = {}
    
    for i in range(num_blocks):
        top_body_name, bottom_body_name = points_body[f'body{i+1}']['bodies'][0].name, points_body[f'body{i+1}']['bodies'][-1].name
        print("Calculating deltas for points between bodies {} and {}".format(top_body_name, bottom_body_name))
        delta_[f'body{i+1}'] = get_deltas(points_body=points_body[f'body{i+1}'], model=model, suspension_body=suspension_body)

    return delta_

def get_linpaths(points_list, model):
    
    
    
    points = points_list['points']

    n = len(points)//4 - 1
    
    linpaths_ = []
    for j in range(4):
        for i in range(n):
            m = 4*i + j
            points_ = (points[m], points[m+4])
            print("Finding linpath between {} and {}".format(points_[0], points_[1]))
            linpath = LinearPathway(points_[0],points_[1])
            linpaths_.append(linpath)
            
    
    bodies = points_list['bodies']
    zipped_bods = []
    for j in range(4):
        for i in range(n):
            m = 4*i + j
            zipped_bods.append((bodies[m], 
                                    bodies[m+4]))
    
    paths_dict = {'paths':linpaths_, 'bodies':zipped_bods}
    
    return paths_dict # linpaths_ # paths_dict

def linpaths_dict(points_body, num_blocks, model):

    linpath_dict = {}
    
    for i in range(num_blocks):
        linpath = get_linpaths(points_body[f'body{i+1}'], model=model)
        linpath_dict[f'body{i+1}'] = linpath


    return linpath_dict
    
def add_wire_info(point_body, body, deltas, linpaths):

    num_wires =  (len(point_body['points']) - 4)//4
    deltas_ = list(deltas[body].values())
    linpaths_ = linpaths[body]
    
    
    
    wire_dict = {}
    wires = {}
    for j in range(4):
        attach_points_ = []
        attach_bodies_mass = []
        for i in range(num_wires):
            m = 4*i + j
            attach_points = (point_body['points'][m], point_body['points'][m+4])
            #print(attach_points)
            attach_bodies = (point_body['bodies'][m].body.mass, point_body['bodies'][m+4].body.mass)
            attach_points_.append(attach_points)
            attach_bodies_mass.append(attach_bodies)

        print(attach_points_)
    
        wire_dict[f'wire{j+1}'] = {'points':attach_points_, 'masses':attach_bodies_mass, 
                                   'deltas':deltas_[num_wires*j:num_wires*j+num_wires],
                                   'linpaths':linpaths_['paths'][num_wires*j:num_wires*j+num_wires], 
                                   'bodies':linpaths_['bodies'][num_wires*j:num_wires*j+num_wires] }
        #print(points_dict)
    
    return wire_dict

def get_wires(num_blocks, point_body, deltas, linpaths):
    
    wires = {}
    for i in range(num_blocks):
        wire_body = add_wire_info(point_body[f'body{i+1}'], body=f'body{i+1}', deltas=deltas, linpaths=linpaths)
        wires[f'body{i+1}'] = wire_body
        
    return wires

def get_force_wire(body, wire, wire_dict, n_wire, k, masses):
    
    g = symbols('g')

    n = n_wire - 1
    path = wire_dict[body][wire]['linpaths'][n]
    points = [path.attachments[0], path.attachments[1]]
    
    bodies = [wire_dict[body][wire]['bodies'][n][0], 
              wire_dict[body][wire]['bodies'][n][1]]

    #print([i.body.name for i in bodies])
    
    
    wire_masses = [i[-1] for i in wire_dict[body][wire]['masses'][:-1]] + [0]
    wire_masses_ = wire_masses[n:]

    #print(wire_masses_)
    
    delta_wire =  wire_dict[body][wire]['deltas'][n_wire-1]
    extension_force = k*delta_wire

    body_gravity = sp.Rational(1, 4) * (np.sum(masses)) * g
    wire_gravity = [i*g for i in wire_masses_]
    
    total_force = [extension_force] + [body_gravity] + wire_gravity
    total_force_ = np.sum(total_force)
    
    return total_force_, points, path, wire_masses_, bodies

def get_forces(n_body, n_wire_body, spring_constants, masses, wire_dict):
    
 
    ks = spring_constants

    wires = {}

    for j in range(4): # 4 is the number of wires per body
        for k in range(n_wire_body): # n_wire_body is the number of sub_wires per body (1 + # of violin mode bodies)
            #print(k+1)
            force = get_force_wire(body=f'body{n_body}', wire=f'wire{j+1}', wire_dict=wire_dict, n_wire=k+1,  k=ks[n_body-1], masses=masses[n_body*5:])[0]
            points = get_force_wire(body=f'body{n_body}', wire=f'wire{j+1}',  wire_dict=wire_dict,  n_wire=k+1,  k=ks[n_body-1],
                                   masses=masses[n_body*5:])[1]
            path = get_force_wire(body=f'body{n_body}', wire=f'wire{j+1}',  wire_dict=wire_dict,  n_wire=k+1,  k=ks[n_body-1],
                                 masses=masses[n_body*5:])[2]
            bodies = get_force_wire(body=f'body{n_body}', wire=f'wire{j+1}',  wire_dict=wire_dict,  n_wire=k+1,  k=ks[n_body-1],
                                   masses=masses[n_body*5:])[-1]

            wires[f'wire{k+1}{j+1}'] = {'force':force, 'points':points, 'path':path, 'bodies':bodies}
        
    return wires


def apply_force_on_body(body_force, suspension_body):

    S = suspension_body
    for i in body_force.keys():
        force = body_force[i]['force']
        path = body_force[i]['path']
        points = body_force[i]['points']
        bodies = body_force[i]['bodies']


        j=0
        for body in bodies:
            body.body.apply_force(path.to_loads(-force, frame = S.global_frame)[j][1], point=points[j])
            j+=1
            #print(body.body.name)
            #print(points[j])

    return 

def apply_gravity(body):
    g = symbols('g')
    body.Fz = -g*body.body.mass
         
    return

    