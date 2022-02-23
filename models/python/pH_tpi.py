# -*- coding: utf-8 -*-
"""
Tellurium oscillation
"""
# Modelling aldolase reaction, pH7.90, DL8

import tellurium as te
import roadrunner
import antimony
import matplotlib.pyplot as plt
# -> reversible
# => irreversible

r = te.loada ('''
model feedback()
  // Reactions:
  vALD: FBP -> GAP + DHAP; 0;
  vGPD: DHAP + NADH -> G3P + NAD; p_GPD1_Vm * (DHAP * NADH - G3P * NAD / p_GPD1_Keq);
  vTPI: DHAP-> GAP; (p_TPI1_Vm / p_TPI1_Kdhap * (DHAP - GAP/p_TPI1_Keq)) / (1 + DHAP/p_TPI1_Kdhap + GAP/p_TPI1_Kglyceral3p);

  // Species initializations:
  FBP = 0; 
  DHAP = 0; 
  GAP = 5.8; 
  G3P = 0; 
  NADH = 0.066; 
  NAD = 0;

  // Variable initialization:
  p_TPI1_Kdhap = 9.5624; 
  p_TPI1_Kglyceral3p = 5.2322; 
  p_TPI1_Vm = 0.0013; 
  p_TPI1_Keq = 0.1092;
  p_GPD1_Vm = 21.7742; 
  p_GPD1_Keq = 3.7037e+04;
  
end''')

result = r.simulate(0, 50, 100)
# r.plot(result)

# selected plot
time = result[:,0]
nadh = result[:,4]
plt.plot(time,nadh)
plt.show()

# print(r.getSBML())