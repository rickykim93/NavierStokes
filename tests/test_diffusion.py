import unittest
import numpy as np
import NavierStokes as ns

class Test_test(unittest.TestCase):

    def test(self):
        # Set the grid dimensions
        lx = 1.
        ly = 1.
        # Set the grid resolution (number of nodes including boundaries)
        nx = 5
        ny = 5
        # The west boundary part, Tw=10 degree
        Tw = 10
        # The east boundary part, Te=100 degree
        Te = 100
        ans=[40.99, 82.45, 108.71, 120.47, 115.17, 29.44, 62.05, 85.38, 98.68, 101.84, 24.44, 50.92, 72.08, 87.03,
             96.35, 22.09, 45.12, 64.98, 81.03, 93.85, 21.09, 42.5, 61.7, 78.25, 92.79]
        test=ns.diffusion(lx, ly, nx, ny, Tw, Te)[1]
        for i in range(5):
            for j in range(5):
                self.assertEqual(np.round(test[i][j], 2), ans[i*5+j])

if __name__ == '__main__':
    unittest.main()
