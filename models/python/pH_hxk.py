# -*- coding: utf-8 -*-
"""
Hexokinase
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
  vGLK: GLCi + ATP -> G6P + ADP; (p_HXK1_Vm * (ATP * GLCi - ( ADP * G6P / p_HXK1_Keq))) / ((p_HXK1_Katp * p_HXK1_Kglc) * (1 + ATP / p_HXK1_Katp + ADP/p_HXK1_Kadp) * (1 + GLCi / p_HXK1_Kglc + G6P / p_HXK1_Kg6p + 0 / p_HXK1_Kt6p));
  vG6PDH: G6P + NADP -> PG6 + NADPH; p_G6PDH_Vm * (NADP * G6P - NADPH * PG6 / p_G6PDH_Keq);

  // Species initializations:
  GLCi = 10;
  G6P = 0;
  PG6 = 0;
  ATP = 1;
  ADP = 0;
  NADP = 1;
  NADPH = 0.0636;

  // Variable initialization:
  p_HXK1_Kadp = 183.8396;
  p_HXK1_Katp = 8.3175e-04;
  p_HXK1_Kg6p = 6.8891e+03;
  p_HXK1_Kglc = 0.0137;
  p_HXK1_Vm = 2.6921e-04;
  p_HXK1_Kt6p = 0.2000;
  p_HXK1_Keq = 7000;
  p_G6PDH_Keq = 16.3000;
  p_G6PDH_Vm = 21.7742;
  
end''')

result = r.simulate(0, 300, 100)
# r.plot(result)

# selected plot
time = result[:,0]
nadh = result[:,7]
plt.plot(time,nadh)
plt.show()

# print(r.getSBML())
