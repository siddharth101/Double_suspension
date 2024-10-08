from utils import (
    add_points, make_wire_bodies, get_points, get_deltas, get_linpaths,
    get_wires, get_forces, deltas_dict, linpaths_dict, add_wire_info,
    get_force_wire, apply_force_on_body, apply_gravity
)
from sympy.physics.mechanics import (
    dynamicsymbols, ReferenceFrame, Point, RigidBody, Particle, System, Force
)
from sympy import zeros, symbols
import sympy as sp
import numpy as np
import sympy.physics.mechanics as me
from IPython.display import display
from suspycious import Model
import suspycious.components as scmp
from pathway import LinearPathway

# Define symbols for distances and weights
d1, d2, n1, n2, w1, w2 = symbols('d1 d2 n1 n2 w1 w2', real=True, positive=True)
d3, n3, w3 = symbols('d3 n3 w3', real=True, positive=True)
d4, n4, w4 = symbols('d4 n4 w4', real=True, positive=True)
d5, n5, w5 = symbols('d5 n5 w5', real=True, positive=True)
d6, n6, w6 = symbols('d6 n6 w6', real=True, positive=True)
d7, n7, w7 = symbols('d7 n7 w7', real=True, positive=True)
d8, d9 = symbols('d8 d9', real=True, positive=True)
l1, l2, l3, l4, l2a = symbols('l1 l2 l3 l4 l_2a', positive=True)


def set_up_bodies():
    """Set up the model and bodies with global positions and orientations."""
    # Create the model and rigid bodies
    model2 = Model()
    S = model2.add(scmp.RigidBody("S"))
    B = model2.add(scmp.RigidBody("B"))
    C = model2.add(scmp.RigidBody("C"))
    D = model2.add(scmp.RigidBody("D"))
    F = model2.add(scmp.RigidBody("F"))

    # Set global positions for the bodies
    S.set_global_position(0, 0, 0)
    B.set_global_position(0, 0, -l1 - d1)
    C.set_global_position(0, 0, -l1 - d1 - d2 - l2 - d3)
    D.set_global_position(0, 0, -l1 - d1 - d2 - l2 - d3 - d4 - l3 - d5)
    F.set_global_position(0, 0, -l1 - d1 - d2 - l2 - d3 - d4 - l3 - d5 - d6 - l4 - d7)

    # Set global orientations for the bodies
    S.set_global_orientation(0, 0, 0)
    B.set_global_orientation(0, 0, 0)
    C.set_global_orientation(0, 0, 0)
    D.set_global_orientation(0, 0, 0)
    F.set_global_orientation(0, 0, 0)

    # Add attachment points to the bodies
    add_points(body=S, point='P', attachment_points=[w1, n1, 0])
    add_points(body=B, point='B', attachment_points=[w1, n1, d1])
    add_points(body=B, point='C', attachment_points=[w2, n2, -d2])
    add_points(body=C, point='D', attachment_points=[w2, n2, d3])
    add_points(body=C, point='E', attachment_points=[w3, n3, -d4])
    add_points(body=D, point='F', attachment_points=[w3, n3, d5])
    add_points(body=D, point='G', attachment_points=[w4, n4, -d6])
    add_points(body=F, point='H', attachment_points=[w4, n4, d7])
    add_points(body=F, point='J', attachment_points=[w5, n5, -d8])

    return S, B, C, D, F, model2


# Set up the bodies and model
S, B, C, D, F, model2 = set_up_bodies()

# Define variables for the number of violin modes
n_violin_modes = 1

# Create body and attachment point lists
bodies_ = [[S, B], [B, C], [C, D], [D, F]]
attach_points_ = [[S.P1, B.B1], [B.C1, C.D1], [C.E1, D.F1], [D.G1, F.H1]]


def get_kane(
    n_violin_modes=n_violin_modes, bodies_=bodies_,
    attach_points_=attach_points_, model=model2, suspension_body=S
):
    """Calculate the dynamics using Kane's method."""
    n_wire_per_wire = n_violin_modes + 1

    points_bodies = get_points(
        n=n_violin_modes, num_blocks=4, model=model,
        bodies=bodies_, attach_points=attach_points_, suspension_body=suspension_body
    )

    deltas = deltas_dict(
        points_body=points_bodies, num_blocks=4, model=model,
        suspension_body=suspension_body
    )

    linpaths = linpaths_dict(
        points_body=points_bodies, num_blocks=4, model=model
    )

    num_blocks = 4
    all_masses = [
        i.body.mass for j in range(num_blocks)
        for i in points_bodies[f'body{j + 1}']['bodies']
    ]

    all_masses_ = list(dict.fromkeys(all_masses))
    wires = get_wires(
        num_blocks=4, point_body=points_bodies,
        deltas=deltas, linpaths=linpaths
    )

    ks = symbols('k1:5')

    body1_force = get_forces(
        n_body=1, n_wire_body=n_wire_per_wire,
        wire_dict=wires, spring_constants=ks, masses=all_masses_
    )
    body2_force = get_forces(
        n_body=2, n_wire_body=n_wire_per_wire,
        wire_dict=wires, spring_constants=ks, masses=all_masses_
    )
    body3_force = get_forces(
        n_body=3, n_wire_body=n_wire_per_wire,
        wire_dict=wires, spring_constants=ks, masses=all_masses_
    )
    body4_force = get_forces(
        n_body=4, n_wire_body=n_wire_per_wire,
        wire_dict=wires, spring_constants=ks, masses=all_masses_
    )

    apply_force_on_body(body1_force, suspension_body=S)
    apply_force_on_body(body2_force, suspension_body=S)
    apply_force_on_body(body3_force, suspension_body=S)
    apply_force_on_body(body4_force, suspension_body=S)

    apply_gravity(model.components['S'])
    apply_gravity(model.components['B'])
    apply_gravity(model.components['C'])
    apply_gravity(model.components['D'])
    apply_gravity(model.components['F']) ### Forgot this earlier

    print("The bodies have been set up, forces applied, now solving for dynamics")
    kane = model.extract_statespace()

    return kane
