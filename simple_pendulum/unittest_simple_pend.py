from sympy import Matrix
from sympy import zeros, symbols
g, m1, k1, l = symbols('g m1 k1 l', positive=True)
A = Matrix([[0,0,1,0],[0,0,0,1],[-g/l, 0,0,0], [0,k1/m1,0,0]])
B = Matrix([[0,0], [0,0], [-g/l, 0], [0, k1/m1]])

import unittest
from simple_pend import simple_pendulum

class test_simple_pendulum(unittest.TestCase):
    def testAmat(self):
        self.assertEqual(simple_pendulum()[0], A)
    def testBmat(self):
        self.assertEqual(simple_pendulum()[1], B)
        
test_simple_pendulum().testAmat()
test_simple_pendulum().testBmat()

if __name__ == "__main__":
    unittest.main()