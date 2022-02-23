# -*- coding: utf-8 -*-
"""
Phosphoglucomutase
"""
# Modelling aldolase reaction, pH7.90, DL4

import tellurium as te
import roadrunner
import antimony
import matplotlib.pyplot as plt
# -> reversible
# => irreversible

r = te.loada ('''
model feedback()
  // Reactions:
  v_PGM: P3G -> P2G; ((p_GPM1_vm / p_GPM1_K3pg) * (P3G - P2G / p_GPM1_Keq)) / (1 + P3G / p_GPM1_K3pg + P2G / p_GPM1_K2pg);
  v_ENO: P2G -> PEP; p_ENO_Vm * (P2G - PEP / p_ENO_Keq);
  v_PYK: PEP + ADP -> PYR + ATP; p_PYK_Vm * (PEP * ADP - PYR * ATP / p_PYK_Keq);
  v_LDH: PYR + NADH -> LAC + NAD; p_LDH1_Vm * (PYR * NADH - LAC * NAD / p_LDH1_Keq);

  // Species initializations:
  P3G = 5;
  P2G = 0;
  PEP = 0;
  PYR = 0;
  LAC = 0;
  ADP = 10;
  ATP = 0;
  NADH = 0.0651;
  NAD = 0;
  
  // Variable initialization:
  p_GPM1_K2pg = 0.0800;
  p_GPM1_K3pg = 1.2031;
  p_GPM1_vm = 0.0023;
  p_GPM1_Keq = 0.1890;
  p_ENO_Vm = 21.7742;
  p_ENO_Keq = 5.1900;
  p_PYK_Vm = 21.7742;
  p_PYK_Keq = 1.0417e+04;
  p_LDH1_Vm = 21.7742;
  p_LDH1_Keq = 8.3333e+03;
  
end''')

result = r.simulate(0, 300, 100)
# r.plot(result)

# selected plot
time = result[:,0]
nadh = result[:,7]
plt.plot(time,nadh)
plt.show()

# print(r.getSBML())

