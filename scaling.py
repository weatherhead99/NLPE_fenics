#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 13:50:43 2016

@author: weatherill
"""

from semicond_quantities import *
#import fipy as fp

from scipy.constants import k, e

class SemiconductorScaling(object):
    def __init__(self, Nmax, Tmin,material=Silicon(),):
        
        self.Ld = DebyeLength(material).calculate(Tmin,Nmax)
        self.ni = IntrinsicConcentration(material).calculate(Tmin)
        self.Vt = k * Tmin / e
        
        self.material = material
#    def get_debye_scaled_spatial_grid(self, Lxyz, debye_scaling):
#        """return a grid with scaling in terms of Debye lengths (natural units)
#            """        
#        if not hasattr(type(Lxyz),'len'):
#            #single value entered, assume 1d
#            meshtp = fp.Grid1D
#            Lxyz = [Lxyz]
#            
#        elif len(Lxyz) == 2:
#            meshtp = fp.Grid2D
#            
#        elif len(Lxyz) == 3:
#            meshtp = fp.Grid3D
#
#        nx = [ _ / (self.Ld * debye_scaling) for _ in Lxyz]
#        dx = [debye_scaling for _ in Lxyz]
#        
#        self.mesh = meshtp(*(dx+ nx))
#
#        return self.mesh

    def scale_V(self,V):
        """scale a voltage (V) to a voltage in Vtherms"""
        return V/ self.Vt
        
    def scale_n(self,n):
        """scale a charge density (cm^-3) to intrinsic concentration per cubic Debye Lengths"""
        return n / self.ni * self.Ld**3
        
    def scale_permittivity(self,eps=None):
        """scale a permittivity (F m^-1) to intrinsic concentration and Debye Length"""            
        if eps is None:
            eps = self.material.CONSTANTS['eps'] * epsilon_0

        eps_um = eps * 1E2 / e #now in electrons cm^-1 Volt^-1
        eps_nat = eps_um / self.ni * self.Vt * self.Ld
        return eps_nat
        
    def scale_distance(self,x):
        """scale a distance (um) to Debye Length"""
        return x / self.Ld

    def scale_mobility(self,mu):
        """scale a mobility in cm**2 / (Vs) to Debye Length and Vtherm"""
        mup = mu / self.Ld**2 * self.Vt
        return mup
        
    
            