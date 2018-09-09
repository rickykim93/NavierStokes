import unittest

from NavierStokes.test import test

class Test_test(unittest.TestCase):

    def test(self):
        self.assertEqual(test(),'test')
    def test2(self):
        self.assertEqual(test(),'t2')

if __name__ == '__main__':
    unittest.main()
