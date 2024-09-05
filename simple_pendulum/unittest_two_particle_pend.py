from sympy import Matrix, eye, zeros, symbols
import unittest
from make_particle_pend import particle_pend

m1, m2, M, g, k1, k2, k3 = symbols('m1, m2, M, g, k1, k2, k3')
L = symbols('L', positive=True)

mat_a1 = Matrix([[-((3*g*(M + m2)/L) + (3*g*(M + m1 + m2)/L))/m1, 0 , 3*g*(M + m2)/(L*m1), 0,0,0],
       [0, -(k1 + k2)/m1, 0, k2/m1, 0,0],
       [3*g*(M + m2)/(L*m2), 0, -(3*g*(M) + 3*g*(M + m2))/(L*m2), 0, 3*M*g/(L*m2), 0],
       [0, k2/m2, 0, -(k2+k3)/m2, 0, k3/m2],
       [0,0, 3*g/L, 0, -3*g/L, 0],
       [0,0,0, k3/M, 0, -k3/M]])

mat_a2 = mat_a1.col_insert(6, zeros(6))

mat_zi = zeros(6).col_insert(6, eye(6))

mat_a = mat_zi.row_insert(6, mat_a2)

mat_a_vals = [mat_a[i,j].simplify() for i in range(12) for j in range(12)]

sys_part_A = particle_pend(n=2)[0]

sys_part_A_vals = [sys_part_A[i,j].simplify() for i in range(12) for j in range(12)]

class test_particle_pendulum(unittest.TestCase):
    def testAmat(self):
        self.assertEqual(sys_part_A_vals, mat_a_vals)
  
        
test_particle_pendulum().testAmat()

if __name__ == "__main__":
    unittest.main()


