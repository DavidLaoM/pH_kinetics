# -*- coding: utf-8 -*-
"""
Pyruvate kinase
"""
# Modelling aldolase reaction, pH7.90, DF4

import tellurium as te
import roadrunner
import antimony
import matplotlib.pyplot as plt
# -> reversible
# => irreversible

r = te.loada ('''
model feedback()
  // Reactions:
  vPYK: ADP + PEP -> ATP + PYR; (((p_PYK1_Vm / (p_PYK1_Kadp * p_PYK1_Kpep) * ADP * PEP) / ((1 + ADP / p_PYK1_Kadp) * (1 + PEP / p_PYK1_Kpep))) * ((PEP / p_PYK1_Kpep + 1) ^ p_PYK1_hill / (p_PYK1_L * ((ATP / p_PYK1_Katp + 1) / (FBP / p_PYK1_Kf16bp + 1)) ^ p_PYK1_hill + (PEP / p_PYK1_Kpep + 1)^p_PYK1_hill)));
  vLDH: NADH + PYR -> NAD + LAC; p_LDH1_Vm * (PYR * NADH - LAC * NAD / p_LDH1_Keq);

  // Species initializations:
  ADP = 10; 
  NADH = 0.0655;
  FBP = 1;
  PEP = 2;
  ATP = 0;
  NAD = 0;
  PYR = 0;
  LAC = 0;

  // Variable initialization:
  p_PYK1_Kadp = 0.2165;
  p_PYK1_Katp = 17.0475;
  p_PYK1_Kf16bp = 0.1892;
  p_PYK1_Kpep = 0.2523;
  p_PYK1_L = 5.6592e+04;
  p_PYK1_hill = 4.2326;
  p_PYK1_Vm = 0.0011;
  p_PYK1_Keq = 1.0417e+04;
  p_LDH1_Vm = 21.7742;
  p_LDH1_Keq = 8.3333e+03;
  
  
end''')

result = r.simulate(0, 600, 100)
# r.plot(result)

# selected plot
time = result[:,0]
nadh = result[:,5]
plt.plot(time,nadh)
plt.show()

# print(r.getSBML())
