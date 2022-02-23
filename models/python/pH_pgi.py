# -*- coding: utf-8 -*-
"""
Phosphoglucoisomerase
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
  v_PGI: G6P -> F6P; (p_PGI1_Vm *(G6P - (F6P / p_PGI1_Keq)))/(1 + G6P / p_PGI1_Kg6p + 1 + F6P / p_PGI1_Kf6p - 1);
  v_PFK: F6P + ATP -> FBP + ADP; p_PFK1_Vm * (F6P * ATP - FBP * ADP / p_PGI1_Keq);
  v_ALD: FBP -> GAP + DHAP; p_FBA1_Vm * (FBP - (GAP * DHAP) / p_FBA1_Keq);
  v_TPI: DHAP -> GAP; p_TPI1_Vm * (DHAP - GAP / p_TPI1_Keq);
  v_GPD: DHAP + NADH -> G3P + NAD; p_GPD1_Vm * (DHAP * NADH - G3P * NAD / p_GPD1_Keq);

  // Species initializations:
  G6P = 5;
  F6P = 0;
  ATP = 1;
  FBP = 0;
  ADP = 0;
  DHAP = 0;
  GAP = 0;
  NADH = 0.0475;
  G3P = 0;
  NAD = 0;
  
  // Variable initialization:
  p_PGI1_Kg6p = 1.2033;
  p_PGI1_Kf6p = 2.9277;
  p_PGI1_Vm = 7.3826e-04;
  p_PGI1_Keq = 0.3660;
  p_PFK1_Vm = 21.7742;
  p_PFK1_Keq = 977;
  p_FBA1_Vm = 21.7742;
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
nadh = result[:,8]
plt.plot(time,nadh)
plt.show()

# print(r.getSBML())