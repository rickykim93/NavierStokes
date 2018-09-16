import unittest
import numpy as np
import NavierStokes.TDMAsolver as TDMA

class Test_test(unittest.TestCase):

    def test(self):
        ##For checking purposes using the example given in class
        aa = np.array([-5., -5., -5., -5.])  # diagonal column below the middle
        bb = np.array([20., 15., 15., 15., 10.])  # middle diagonal column value
        cc = np.array([-5., -5., -5., -5.])  # diagonal column above the middle
        dd = np.array([1100., 100., 100., 100., 100.])  # the right column (b)
        self.assertEqual(np.round(TDMA.TDMAsolver(aa, bb, cc, dd),2)[0],64.23)
        self.assertEqual(np.round(TDMA.TDMAsolver(aa, bb, cc, dd),2)[1],36.91)
        self.assertEqual(np.round(TDMA.TDMAsolver(aa, bb, cc, dd),2)[2],26.50)
        self.assertEqual(np.round(TDMA.TDMAsolver(aa, bb, cc, dd),2)[3],22.60)
        self.assertEqual(np.round(TDMA.TDMAsolver(aa, bb, cc, dd),2)[4],21.30)

if __name__ == '__main__':
    unittest.main()
