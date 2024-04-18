# %%
import sympy as sp
import sympy.physics.mechanics as me
from suspycious import Model
import suspycious.components as scmp
import suspycious.constants as const
me.init_vprinting()
# %%
L, k = sp.symbols("L, k", positive=True)
d = sp.symbols("d", real=True)
# %%
model = Model()
A = model.add(scmp.RigidBody("A"))
S = model.add(scmp.RigidBody("S"))
A.set_global_position(0, 0, -L)
A.set_global_orientation(0, 0, 0)
S.set_global_position(0, 0, 0)
S.set_global_orientation(0, 0, 0)
# %%
P = A.add_fixed_point("P", 0, 0, d)
# %%
dx = S.com.pos_from(P).express(model.frame)
dx0 = model.op_subs(dx)
# dx0 = model.op_subs(P.pos_from(S.com).express(model.frame))
# dx = P.pos_from(S.com)
# ddx = dx - dx0
# ddx = ddx.express(model.frame)
# %%
# A.body.apply_force(k * P.pos_from(S.com).express(model.frame))
A.body.apply_force(force=k * (dx - dx0) + A.M * const.g * model.frame.z, point=P)
# A.body.apply_force(force=k * (dx - dx0) + A.M * const.g * dx.normalize(), point=P)
# A.body.apply_force(force=k * (dx - dx0))
model.apply_gravity()
# A.body.apply_force(-)
# %%
# d_from_origin = P.pos_from(model.origin)
# display(d_from_origin)
# d_from_origin = d_from_origin.express(model.frame).subs(A.op_point)
# display(d_from_origin)
# %%
ss = model.extract_statespace(include_velocities=True)
# %%
eomsA = ss.C @ ss.A @ ss._state_vector
eomsB = ss.C @ ss.B @ ss._input_vector
# %%
display(eomsA[ss.outputs[A.ux], 0])
display(eomsB[ss.outputs[A.ux], 0])
# %%
display(eomsA[ss.outputs[A.wy], 0])
display(eomsB[ss.outputs[A.wy], 0])
# %%
display(eomsA[ss.outputs[A.wx], 0])
display(eomsB[ss.outputs[A.wx], 0])
# %%
display(eomsA[ss.outputs[A.uz], 0])
display(eomsB[ss.outputs[A.uz], 0])
# %%
