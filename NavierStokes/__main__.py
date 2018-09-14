import argparse
import NavierStokes as NS
from NavierStokes.TDMAsolver import main as TDMA

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-td","--tdma", help="TDMA solver", action="store_true")
    parser.add_argument("-t", help="just a test")
    parser.add_argument("-v", "--version", help="version", action="store_true")
    a = parser.parse_args()

    if a.version:
        return NS.__version__
    elif a.tdma:
        TDMA()
