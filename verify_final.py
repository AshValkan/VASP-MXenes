#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify the NEB final-state POSCAR: report per-Li displacement vs the relaxed initial state.
Usage: python3 verify_final.py <initial_CONTCAR> <final_POffset POSCAR>"""
import sys, numpy as np

def load(path):
    L = open(path).read().splitlines()
    s = float(L[1]); A = np.array([[float(x) for x in L[2+i].split()[:3]] for i in range(3)])*s
    sp = L[5].split(); ct = [int(x) for x in L[6].split()]
    mi = next(i for i in range(7,12) if L[i].strip() and L[i].strip()[0].upper() in 'DCS')
    if L[mi].strip()[0].upper()=='S': mi += 1
    assert L[mi].strip()[0].upper()=='D', 'expect Direct coords'
    n = sum(ct)
    D = np.array([[float(x) for x in L[mi+1+i].split()[:3]] for i in range(n)])
    return A, sp, ct, D@A

A, sp, ct, C0 = load(sys.argv[1])
_, _, _, C1 = load(sys.argv[2])
assert sp[-1] == 'Li', f'last species expected Li, got {sp[-1]}'
n = sum(ct); i0 = n - ct[-1]
for j in range(ct[-1]):
    d = np.linalg.norm(C1[i0+j]-C0[i0+j])
    print(f'Li#{i0+j}: init {C0[i0+j].round(3)} -> final {C1[i0+j].round(3)}  |d|={d:.3f} A')
