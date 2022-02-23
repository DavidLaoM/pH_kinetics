# -*- coding: utf-8 -*-
"""
Pyruvate decarboxylase
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
  vPDC: PYR -> CO2 + AcAld; (p_PDC1_vm * (PYR / p_PDC1_Kpyr) ^ p_PDC1_hill) / (1 + (PYR / p_PDC1_Kpyr) ^ p_PDC1_hill);
  vADH: AcAld + NADH -> ETOH + NAD; p_ADH_Vm * (AcAld * NADH - ETOH * NAD / p_ADH_Keq);

  // Species initializations:
  PYR = 50;
  CO2 = 0;
  AcAld = 0;
  ETOH = 0;
  NADH = 0.0823;
  NAD = 0;

  // Variable initialization:
  p_PDC1_Kpyr = 8.2185;
  p_PDC1_hill = 2.0231;
  p_PDC1_vm = 7.7593e-04;
  p_ADH_Vm = 21.7742;
  p_ADH_Keq = 769.2308;
  
end''')

result = r.simulate(0, 300, 100)
# r.plot(result)

# selected plot
time = result[:,0]
nadh = result[:,4]
plt.plot(time,nadh)
plt.show()

# print(r.getSBML())
