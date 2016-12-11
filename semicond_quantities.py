#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 13:10:19 2016

@author: weatherill
"""
from __future__ import division

from scipy.constants import e, k, m_e, pi, h, epsilon_0
import math

class Semiconductor(object):
    def abs_permittivity(self):
        return self.CONSTANTS['eps'] * epsilon_0
        

class Silicon(Semiconductor):
    CONSTANTS = {"Eg" : 1.12 # ev
                 , "m*eDOS" : 1.08 #relative to m_e
                 , "m*hDOS" : 0.81
                 , 'eps' : 11.68
                 }
    

class SiO2(Semiconductor):
    CONSTANTS = {'eps' : 3.9}
                 
class SemiconductorQuantity(object):
    def __init__(self,material=Silicon()):
        self.material = material
        


class IntrinsicConcentration(SemiconductorQuantity):
    def calculate(self,T):
        """ return ni (in cm^-3) for temperature T"""
        #DOS in the conduction band
        Nc = 2 * ( 2 * self.material.CONSTANTS['m*eDOS'] * m_e* pi * k * T /h**2)**1.5
        #DOS in the valend band
        Nv = 2 * ( 2 * self.material.CONSTANTS['m*hDOS'] *m_e* pi * k * T / h**2)**1.5

        Eg = self.material.CONSTANTS['Eg'] * e
        ni = math.sqrt(Nc * Nv) * math.exp(-Eg / (2* k *T ) ) # m^-3
        
        return ni * 1E-6 #cm^-3


class DebyeLength(SemiconductorQuantity):
    def calculate(self,T,N):
        """return Debye Length (in um) for temperature T, doping N"""
        
        Np = N*1E6 #convert from cm^-3 to m^-3
        
        Ld  = math.sqrt( self.material.CONSTANTS['eps'] * epsilon_0 * k * T / (e**2 * Np))
        return Ld * 1E6 #um

class ElectronMobilityCaugheyThomas(SemiconductorQuantity):
    P = {"b0i" : 1.109,
         "b1i" : 0.66,
         "Tbi" : 300,
         "v0i" : 2.4E7,
         "thvi": 0.8,
         "Tvi" : 600,
         "mui1" : 55.24,
         "mui2" : 1429.23,
         "mui3" : -2.3,
         "mui4" : -3.8,
         "mui5" : 0.73,
         "Tmui" : 300,
         "Nmui" : 1.073E17}
    def calculate(self,T,N,Efield):
        betai = self.P['b0i'] * (T / self.P['Tbi']) ** self.P['b1i']
        vsati = self.P['v0i'] / ( 1 + self.P['thvi'] * math.exp(T / self.P['Tvi']))
        
        Tui = T / self.P['Tmui']
        
        mui0 = self.P['mui1'] + (self.P['mui2']*Tui**self.P['mui3'] - self.P['mui1'] ) \
                / ( 1 + Tui**self.P['mui4'] * (N/self.P['Nmui'])**self.P['mui5'])

        mu = mui0 / ( 1 + (mui0*Efield / vsati)**betai  )**(1./betai)
        
        return mu
        
class HoleMobilityCaugheyThomas(ElectronMobilityCaugheyThomas):
    def __init__(self):
        self.P = self.P.copy()
        self.P['b0i'] = 1.213
        self.P['b1i'] = 0.17
        self.P['mui1'] = 49.7
        self.P['mui2'] = 479.37
        self.P['mui3'] = -2.2
        self.P['mui4'] = -3.7
        self.P['mui5'] = 0.70
        self.P['Nmui'] = 1.606E17


