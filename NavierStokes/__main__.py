import argparse
from NavierStokes.test import test
from NavierStokes.TDMAsolver import TDMAsolver

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-t", help="just a test")
    a = parser.parse_args()

    if a.t:
        test()