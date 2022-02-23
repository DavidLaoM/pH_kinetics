# -*- coding: utf-8 -*-
"""
Aldolase
"""
# Modelling aldolase reaction, pH7.90, DF1

import tellurium as te
import roadrunner
import antimony
import matplotlib.pyplot as plt
# -> reversible
# => irreversible

r = te.loada ('''
model feedback()
  // Reactions:
  vALD: FBP -> GAP + DHAP; p_FBA1_Vm * (FBP - (GAP * DHAP) / p_FBA1_Keq) / (p_FBA1_Kf16bp * (FBP / p_FBA1_Kf16bp + (1 + GAP / p_FBA1_Kglyceral3p) * (1 + DHAP / p_FBA1_Kdhap)));
  vTPI: DHAP-> GAP; p_TPI1_Vm * (DHAP - GAP / p_TPI1_Keq);
  vGPD: DHAP + NADH -> G3P + NAD; p_GPD1_Vm * (DHAP * NADH - G3P * NAD / p_GPD1_Keq);

  // Species initializations:
  FBP = 2; 
  DHAP = 0; 
  GAP = 0; 
  G3P = 0; 
  NADH = 0.0945; 
  NAD = 0;

  // Variable initialization:
  p_FBA1_Kf16bp = 0.4303; 
  p_FBA1_Kglyceral3p = 632.4588; 
  p_FBA1_Kdhap = 0.0074; 
  p_FBA1_Vm = 7.2458e-04; 
  p_FBA1_Keq = 6.5000e-04;
  p_TPI1_Vm = 21.7742; 
  p_TPI1_Keq = 0.1065;
  p_GPD1_Vm = 21.7742; 
  p_GPD1_Keq = 5000;
  
end''')

result = r.simulate(0, 400, 100)
# r.plot(result)

# selected plot
time = result[:,0]
nadh = result[:,4]
plt.plot(time,nadh)
plt.show()

# print(r.getSBML())
