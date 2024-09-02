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


def get_linpaths(body1, body2, points_1='top', points_2='top', suspension_body=None):
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
    linpaths = {
        'paths': [LinearPathway(i, j) for i, j in zip(body1_points, body2_points)],
        'bodies': [body1.body, body2.body]
    }

    return linpaths


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