# -*- coding: utf-8 -*-
"""
FILE TDMAsolver.py
Solves TDMA

Name: Kyu Mok Kim
Student Number: 998745381

MIE1210 Project 1

@author: kimkyu4

This code is designed to solve the system of equations denoted as: Ax=b,
using TDMA method using a,b,c,d arrays.
Refer to my report for detailed explanation.
"""

import numpy as np
import xml.etree.ElementTree as ET

# TDMA Solver
def TDMAsolver(aa, bb, cc, dd):
    a, b, c, d = map(np.array, (aa, bb, cc, dd))  # preventing overwriting original
    n = len(b)  # number of rows/equations

    for i in range(1, n):
        e = a[i - 1] / b[i - 1]
        b[i] = b[i] - e * c[i - 1]
        d[i] = d[i] - e * d[i - 1]

    x = b
    x[-1] = d[-1] / b[-1]

    for j in range(n - 2, -1, -1):
        x[j] = (d[j] - c[j] * x[j + 1]) / b[j]

    return x

def main():
    #TODO fix config.xml reading on module
    tree = ET.parse('config.xml')
    root=tree.getroot()
    for child in root:
        if child.tag=='TDMAsolver':
            return TDMAsolver([float(s) for s in child.find('aa').text.split(',')],
                              [float(s) for s in child.find('bb').text.split(',')],
                              [float(s) for s in child.find('cc').text.split(',')],
                              [float(s) for s in child.find('dd').text.split(',')])

if __name__ == '__main__':
    print (main())