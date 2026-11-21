1 # !/ usr / bin / env python3
2 " " "
3 eht_shadow_check . py
4 Estimate of the consciousness - point correction to the M87 * shadow diameter .
5 " " "
6 import numpy as np
7
8 M = 6.5 e9 * 1.989 e30
9 G = 6.674 e -11
10 c = 2.998 e8
11 lP = 1.616 e -35
12 R_s = 2* G * M / c **2
13 delta = lP / R_s
14 print ( f " M87 * Schwarzschild radius R_s = { R_s :.3 e } m " )
15 print ( f " Fractional correction delta theta / theta ~ l_P / R_s = { delta :.3 e } " )
10
16 print ( " Current EHT precision is far too coarse to detect this correction . " )
