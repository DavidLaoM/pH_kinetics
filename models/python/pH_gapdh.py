# -*- coding: utf-8 -*-
"""
Tellurium oscillation
"""
# Modelling aldolase reaction, pH7.81, D1

import tellurium as te
import roadrunner
import antimony
import matplotlib.pyplot as plt
# -> reversible
# => irreversible

r = te.loada ('''
model feedback()
  // Reactions:
  vPGK: BPG + ADP -> P3G + ATP; p_PGK_Vm * (BPG * ADP - P3G * ATP / p_PGK_Keq);
  vGAPDHfwd: NAD + GAP + PHOS -> BPG + NADH; (p_TDH1_Vmf * (GAP * NAD - BPG * NADH / p_GAPDH_Keq) / (p_TDH1_Kgap * p_TDH1_Knad)) / ((1 + NAD / p_TDH1_Knad + NADH / p_TDH1_Knadh) * (1 + BPG / p_TDH1_Kbpg + GAP / p_TDH1_Kgap));
  
  // Species initializations:
  P3G = 0; 
  ATP = 0; 
  BPG = 0; 
  ADP = 10; 
  NAD = 1; 
  GAP = 5.8; 
  PHOS = 0; 
  NADH = 0;

  // Variable initialization:
  p_TDH1_Kgap = 0.9678; 
  p_TDH1_Kbpg = 0.0559; 
  p_TDH1_Knad = 0.8771; 
  p_TDH1_Knadh = 0.0097; 
  p_TDH1_Vmf = 0.0041; 
  p_GAPDH_Keq = 0.0123;
  p_PGK_Vm = 21.7742; 
  p_PGK_Keq = 1.9231e+03;
  
end''')

result = r.simulate(0, 300, 100)
# r.plot(result)

# selected plot
time = result[:,0]
nadh = result[:,8]
plt.plot(time,nadh)
plt.show()

# print(r.getSBML())