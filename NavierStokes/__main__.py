import argparse
import NavierStokes as ns
from subprocess import call
import os

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-td","--tdma", help="TDMA solver", action="store_true")
    parser.add_argument("-e","--edit", help="edit config file", action="store_true")
    parser.add_argument("-v", "--version", help="version", action="store_true")
    a = parser.parse_args()

    if a.version:
        return ns.__version__
    elif a.tdma:
        ns.cmd_TDMA()
    elif a.edit:
        call('nano {}'.format(os.path.join(os.path.dirname(os.path.realpath(__file__)),'config.xml')), shell=True)