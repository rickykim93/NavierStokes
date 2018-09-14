import argparse
from NavierStokes.test import test
from NavierStokes.TDMAsolver import TDMAsolver

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-td","--tdma", help="TDMA solver", action="store_true")
    parser.add_argument("-t", help="just a test")
    a = parser.parse_args()

    if a.t:
        test()
    elif a.tdma:
        #TODO: make this work
        TDMAsolver()
