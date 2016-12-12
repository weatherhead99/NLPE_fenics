#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 21:43:59 2016

@author: danw
"""


import fenics as fc
from scaling import SemiconductorScaling

Nd = 1E16 #cm^-3
T=270

xmax = 5
xmin = -5

debyefact = 0.3
scaler = SemiconductorScaling(Nd,T)
perm= scaler.scale_permittivity()

#TODO: incorporate in scaling in a mesh-independant way
num_LDs =  abs(xmax - xmin) / scaler.Ld
num_meshpts = int( num_LDs / debyefact)

mesh = fc.IntervalMesh(num_meshpts,scaler.scale_distance(-5),scaler.scale_distance(5))

V = fc.FunctionSpace(mesh,'P',1)


def boundary(x, on_boundary):
    return on_boundary
    
    
bc = fc.DirichletBC(V,fc.Constant(0.),boundary)

potential = fc.Function(V)
t = fc.TestFunction(V)

doping = fc.Expression(' (x[0] < 0) ? N_D : -N_A ', N_D=scaler.scale_n(Nd), N_A=scaler.scale_n(Nd),degree=2 )

n = fc.Expression('exp(pot) * exp(-phi_n)',pot=potential,phi_n=0,degree=2)
p = fc.Expression('exp(-pot) * exp(phi_p)',pot=potential,phi_p=0,degree=2)

init_pot_guess = fc.Expression('asinh(D/2.)',D=doping,degree=2)
potential.assign(init_pot_guess)

dP = fc.TrialFunction(V)

#calculate 1 potential delta
dp_bilinear = perm * fc.dot(fc.grad(dP),fc.grad(t))*fc.dx - dP * t * (n + p)*fc.dx
dp_linear = (n - p - doping - (n+p)*potential)*t*fc.dx


dPP = fc.Function(V)
fc.solve(dp_bilinear == dp_linear, dPP,bc)

