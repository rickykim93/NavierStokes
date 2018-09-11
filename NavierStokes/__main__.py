import argparse
from NavierStokes.test import test
from NavierStokes.TDMAsolver import TDMAsolver

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-aa", help="diagonal column below the middle")
    parser.add_argument("-bb", help="middle diagonal column value")
    parser.add_argument("-cc", help="diagonal column above the middle")
    parser.add_argument("-dd", help="the right column (b)")
    parser.add_argument("-td","--tdma", help="TDMA solver", action="store_true")
    parser.add_argument("-t", help="just a test")
    a = parser.parse_args()

    if a.t:
        test()
    elif a.tdma:
        TDMAsolver(a.aa,a.bb,a.cc,a.dd)
