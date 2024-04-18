# %%
import numpy as np
import sympy as sp
import sympy.physics.mechanics as me

from suspycious import Model
import suspycious.components as scmp

me.init_vprinting()
# %%
kx, ky, kz = sp.symbols("k_x, k_y, k_z", positive=True)
Kx, Ky, Kz = sp.symbols("kappa_x, kappa_y, kappa_z", positive=True)
# %%
model = Model()
A = model.add(scmp.RigidBody("A"))
# %%
# A.body.clear_loads()
A.body.apply_force(-kx * A.x_cm * A.frame.x)
A.body.apply_force(-ky * A.y_cm * A.frame.y)
A.body.apply_force(-kz * A.z_cm * A.frame.z)
A.body.apply_torque(-Kx * A.tx_cm * A.frame.x)
A.body.apply_torque(-Ky * A.ty_cm * A.frame.y)
A.body.apply_torque(-Kz * A.tz_cm * A.frame.z)
# %%
ss = model.extract_statespace(include_velocities=True)
display(ss.E)
display(ss.A)
display(ss.B)
display(ss._sol[-1])
# %%
eomsA = ss.C @ ss.A @ ss._state_vector
eomsB = ss.C @ ss.B @ ss._input_vector
# %%
display(eomsA[ss.outputs[A.ux], 0])
display(eomsB[ss.outputs[A.ux], 0])
# %%
display(eomsA[ss.outputs[A.wx], 0])
display(eomsB[ss.outputs[A.wx], 0])
# %%
