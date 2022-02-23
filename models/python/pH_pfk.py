# -*- coding: utf-8 -*-
"""
Phosphofructokinase
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
  v_PFK: F6P + ATP -> FBP + ADP; (p_PFK_Vm * p_PFK_gR * (F6P / p_PFK_Kf6p) * (ATP / p_PFK_Katp) * (1 + (F6P / p_PFK_Kf6p) + (ATP / p_PFK_Katp) + p_PFK_gR * ((F6P / p_PFK_Kf6p) * (ATP / p_PFK_Katp)))) / ((1 + F6P / p_PFK_Kf6p + ATP / p_PFK_Katp + (p_PFK_gR * (F6P / p_PFK_Kf6p) * (ATP / p_PFK_Katp))) ^ 2 + p_PFK_L * ((1 + p_PFK_Ciatp * (ATP / p_PFK_Kiatp)) / (1 + ATP / p_PFK_Kiatp)) ^ 2 * ((1 + p_PFK_Camp * (AMP / p_PFK_Kamp)) / (1 + AMP / p_PFK_Kamp)) ^ 2 * ((1 + ((p_PFK_Cf26bp * F26BP) / p_PFK_Kf26bp) + ((p_PFK_Cf16bp * FBP) / p_PFK_Kf16bp)) / (1 + (F26BP / p_PFK_Kf26bp) + (FBP / p_PFK_Kf16bp))) ^ 2 * (1 + p_PFK_Catp * (ATP / p_PFK_Katp)) ^ 2);
  v_ALD: FBP -> DHAP + GAP; p_FBA1_Vm * (FBP - (GAP * DHAP) / p_FBA1_Keq);
  v_GPD: DHAP + NADH -> G3P + NAD; p_GPD1_Vm * (DHAP * NADH - G3P * NAD / p_GPD1_Keq);
  v_TPI: DHAP -> GAP; p_TPI1_Vm * (DHAP - GAP / p_TPI1_Keq);

  // Species initializations:
  FBP = 0;
  DHAP = 0;
  GAP = 0;
  G3P = 0;
  NADH = 0.0899;
  NAD = 0;
  F6P = 10;
  ATP = 0.5;
  ADP = 0;

  // Variable initialization:
  p_PFK_gR =20.0132;
  p_PFK_Kf6p =0.0702;
  p_PFK_Katp =0.1607;
  p_PFK_L =0.6246;
  p_PFK_Ciatp =51.2436;
  p_PFK_Kiatp =0.9815;
  p_PFK_Camp =0.0845;
  p_PFK_Kamp =0.0950;
  p_PFK_Cf26bp =0.0155;
  p_PFK_Kf26bp =0.0013;
  p_PFK_Cf16bp =0.3604;
  p_PFK_Kf16bp =0.0949;
  p_PFK_Catp =15.0631;
  p_PFK_F26BP =0.001;
  p_PFK_Vm =4.2029e-05;
  p_PFK_Keq = 977;
  p_FBA1_Keq = 6.5000e-04;  p_FBA1_Vm = 21.7742;
  p_GPD1_Keq = 5000;  p_GPD1_Vm = 21.7742;
  p_TPI1_Keq = 0.1066;  p_TPI1_Vm = 21.7742; 
  F26BP = 0.001;
  AMP = 0;
  
end''')

result = r.simulate(0, 1000, 100)
# r.plot(result)

# selected plot
time = result[:,0]
nadh = result[:,7]
plt.plot(time,nadh)
plt.show()

# print(r.getSBML())
