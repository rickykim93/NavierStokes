import argparse
from NavierStokes.test import test

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("t", help="just a test")
    a = parser.parse_args()

    test()