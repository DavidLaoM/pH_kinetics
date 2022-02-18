% pH values tested
% pH_vals     = [6.19 6.26 6.41 6.60 6.81 7.06 7.29 7.51 7.68 7.81];
pH_vals     = [6.19 6.26 6.32 6.41 6.60 6.81 7.06 7.29 7.51 7.68 7.81 7.90];
% keq manually taken from the equilibrator. All concentrations in [mM]
% units.
kEq.HXK     = [2.74E2 3.03E2 3.31E2 3.79E2 5.16E2 7.45E2 1.2E3 1.9E3 3E3 4.3E3 5.7E3 7.0E3];  %dir+ EC 2.7.1.1 or EC 2.7.1.2
kEq.PGI     = [3.6E-1 3.6E-1 3.6E-1 3.6E-1 3.6E-1 3.61E-1 3.61E-1 3.62E-1 3.63E-1 3.64E-1 3.65E-1 3.65E-1];  %dir+
kEq.PFK     = [2.56E1 2.96E1 3.34E1 4.03E1 5.95E1 9.15E1 1.54E2 2.5E2 4.05E2 5.91E2 7.95E2 9.77E2];  %dir+, E.C.num.: EC 2.7.1.56
kEq.FBA     = [1.2E-3 1.1E-3 1.0E-3 9.7E-4 8.6E-4 7.9E-4 7.3E-4 7E-4 6.8E-4 6.7E-4 6.6E-4 6.5E-4];  %dir+ [EC 4.1.2.13]
kEq.TPI     = [1/(8.07) 1/(8.21) 1/(8.30) 1/(8.46) 1/(8.74) 1/(8.97) 1/(9.16) 1/(9.26) 1/(9.33) 1/(9.36) 1/(9.38) 1/(9.38)];  %dir-
kEq.GAPDH   = [1.7E-3, 2.3E-3, 2.9E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1,1.78E-1, 2.45E-1, 3.05E-1];  %dir+ [EC 1.2.1.13] [EC 1.2.1.59]
kEq.PGK     = [1/(7.4E-4), 1/(7.2E-4), 1/(7.1E-4), 1/(6.9E-4), 1/(6.5E-4), 1/(6.1E-4), 1/(5.7E-4), 1/(5.5E-4), 1/(5.3E-4), 1/(5.2E-4), 1/(5.2E-4), 1/(5.1E-4)];  %dir- [EC 2.7.2.3]
kEq.ENO     = [5.22 5.22 5.22 5.21 5.21 5.2 5.2 5.2 5.19 5.19 5.19 5.19];  %dir ?
kEq.PGM     = [1/6.62 1/6.44 1/6.29 1/6.12 1/5.83 1/5.62 1/5.46 1/5.38 1/5.33 1/5.3 1/5.29 1/5.29]; %dir-
kEq.PYK     = [1/(3.1E-6) 1/(3.5E-6) 1/(3.9E-6) 1/(4.5E-6) 1/(6.5E-6) 1/(9.7E-6) 1/(1.6E-5) 1/(2.5E-5) 1/(4.1E-5) 1/(5.9E-5) 1/(7.9E-5) 1/(9.6E-5)];  %dir-
kEq.PDC     = [1E4 8.9E3 7.8E3 6.3E3 4.1E3 2.5E3 1.4E3 8.31E2 5.01E2 3.39E2 2.51E2 2.04E2];  %dir+ [EC 4.1.1.1]
% For GPD: 1.1.1.8 (EC 1.1.1.8)
kEq.GPD     = [1/(2.9E-6) 1/(3.5E-6) 1/(4.2E-6) 1/(5.3E-6) 1/(8.7E-6) 1/(1.5E-5) 1/(2.7E-5) 1/(4.7E-5) 1/(7.9E-5) 1/(1.2E-4) 1/(1.6E-4) 1/(2.0E-4) ]; %dir-
% For LDH: 1.1.1.27 (EC 1.1.1.27)
kEq.LDH     = [1/(2.3E-6)	1/(2.7E-6) 1/(3.2E-6) 1/(3.9E-6) 1/(6.1E-6)	1/(9.9E-6) 1/(1.8E-5) 1/(3.0E-5) 1/(5.0E-5) 1/(7.3E-5) 1/(9.9E-5) 1/(1.2E-4)];
% For ADH: (EC 1.1.1.1)
full_pH_vals     = [6.19 6.26 6.32 6.41 6.60 6.81 7.06 7.29 7.51 7.68 7.81 7.90];
kEq.ADH     = [1/(2.4E-5) 1/(2.8E-5) 1/(3.3E-5) 1/(4.0E-5) 1/(6.3E-5) 1/(1.0E-4) 1/(1.8E-4) 1/(3.1E-4) 1/(5.1E-4) 1/(7.6E-4) 1/(1.0E-3) 1/(1.3E-3)];

figure(1)

subplot(3,4,1) % HXK
semilogy(pH_vals,kEq.HXK,'.-')
title('HXK')
xlabel('pH')
ylabel('k_{EQ}'), set(gca,'xtick',[])

subplot(3,4,2) % PGI
semilogy(pH_vals,kEq.PGI,'.-')
set(gca,'Color',[0.9 0.9 0.9])
title('PGI')
xlabel('pH')
ylabel('k_{EQ}'), set(gca,'xtick',[])

subplot(3,4,3) % PFK
semilogy(pH_vals,kEq.PFK,'.-')
title('PFK')
xlabel('pH')
ylabel('k_{EQ}'), set(gca,'xtick',[])

subplot(3,4,4) % FBA
semilogy(pH_vals,kEq.FBA,'.-')
title('FBA')
xlabel('pH')
ylabel('k_{EQ}'), set(gca,'xtick',[])

subplot(3,4,5) % TPI
semilogy(pH_vals,kEq.TPI,'.-')
set(gca,'Color',[0.9 0.9 0.9])
title('TPI')
xlabel('pH')
ylabel('k_{EQ}'), set(gca,'xtick',[])

subplot(3,4,6) % GAPDH
semilogy(pH_vals,kEq.GAPDH,'.-')
title('GAPDH')
xlabel('pH')
ylabel('k_{EQ}'), set(gca,'xtick',[])

subplot(3,4,7) % PGK
semilogy(pH_vals,kEq.PGK,'.-')
title('PGK')
xlabel('pH')
ylabel('k_{EQ}'), set(gca,'xtick',[])

subplot(3,4,8) % ENO
semilogy(pH_vals,kEq.ENO,'.-')
set(gca,'Color',[0.9 0.9 0.9])
title('ENO')
xlabel('pH')
ylabel('k_{EQ}'), set(gca,'xtick',[])

subplot(3,4,9) % PGM
semilogy(pH_vals,kEq.PGM,'.-')
set(gca,'Color',[0.9 0.9 0.9])
title('PGM')
xlabel('pH')
ylabel('k_{EQ}'), set(gca,'xtick',[])

subplot(3,4,10) % PYK
semilogy(pH_vals,kEq.PYK,'.-')
title('PYK')
xlabel('pH')
ylabel('k_{EQ}'), set(gca,'xtick',[])

subplot(3,4,11) % PDC
semilogy(pH_vals,kEq.PDC,'.-')
title('PDC')
xlabel('pH')
ylabel('k_{EQ}'), set(gca,'xtick',[])

suptitle('k_{EQ} glycolysis enzymes vs pH (Source: eQuilibrator)')

set(gcf,'color','w');

%%
figure(2)
semilogy(pH_vals,kEq.HXK,'.-','MarkerSize',20,'MarkerFaceColor','Blue','Color','blue','LineWidth',2)
set(gca,'fontsize',18)
title('HXK','fontsize', 28)
xlabel('pH','fontsize', 28)
ylabel('k_{EQ}','fontsize', 28)
% plot

save('pHchange_kEq.mat','pH_vals','kEq')




