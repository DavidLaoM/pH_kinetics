# -*- coding: utf-8 -*-
"""
Enolase
"""
# Modelling enolase reaction, pH7.90, DF1

import tellurium as te
import roadrunner
import antimony
import matplotlib.pyplot as plt
# -> reversible
# => irreversible

r = te.loada ('''
model feedback()
  // Reactions:
  vENO: P2G -> PEP; ((p_ENO1_vm / p_ENO1_K2pg) * (P2G - PEP / p_ENO1_Keq)) / (1 + P2G / p_ENO1_K2pg + PEP / p_ENO1_Kpep);
  
  // Species initializations:
  P2G = 6; 
  PEP = 0.4975;

  // Variable initialization:
  p_ENO1_K2pg = 6.7187e-04; 
  p_ENO1_Kpep = 64.5850; 
  p_ENO1_vm = 0.0011; 
  p_ENO1_Keq = 5.19;
  
end''')

result = r.simulate(0, 600, 100)
# r.plot(result)

# selected plot
time = result[:,0]
p2g = result[:,2]
plt.plot(time,p2g)
plt.show()

# print(r.getSBML())

