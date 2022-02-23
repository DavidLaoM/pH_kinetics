# -*- coding: utf-8 -*-
"""
Glyceraldehyde dehydrogenase (reverse direction)
"""
# Modelling aldolase reaction, pH7.81, DF1

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
  vGAPDHrev: BPG + NADH -> NAD + GAP + PHOS; (p_TDH1_Vmr * (BPG * NADH - GAP * NAD * p_GAPDH_Keq) / (p_TDH1_Kbpg * p_TDH1_Knadh)) / ((1 + NAD / p_TDH1_Knad + NADH / p_TDH1_Knadh) * (1 + BPG / p_TDH1_Kbpg + GAP / p_TDH1_Kgap));

  // Species initializations:
  P3G = 5;
  ATP = 1;
  BPG = 0;
  ADP = 0;
  NAD = 0;
  GAP = 0;
  PHOS = 500;
  NADH = 0.0875;

  // Variable initialization:
  p_TDH1_Kgap = 0.0222;
  p_TDH1_Kbpg = 0.0143;
  p_TDH1_Knad = 0.0071;
  p_TDH1_Knadh = 5.8235e-04;
  p_TDH1_Vmr = 0.0015;
  p_GAPDH_Keq = 0.0123;
  p_PGK_Vm = 21.7742;
  p_PGK_Keq = 1.9231e+03;
  
end''')

result = r.simulate(0, 300, 100)
# r.plot(result)

# selected plot
time = result[:,0]
nadh = result[:,5]
plt.plot(time,nadh)
plt.show()

# print(r.getSBML())

