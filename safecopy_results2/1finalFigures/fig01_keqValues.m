% % FIG01_KEQVALUES
% In this figure the pH-dependent keq for gapdh and pgk are calculated and 
% displayed.

% pH values
pH_vals =   [6.1900;
            6.2600;
            6.4100;
            6.6000;
            6.8100;
            7.0600;
            7.2900;
            7.5100;
            7.6800;   
            7.8100];
% keq.gapdh.reference
keq_gapdh_ref = 0.029 * 0.05; %@pH=7, added 50 mM PHOS
% keq.gapdh.calculated
keq_gapdh = keq_gapdh_ref .* 10 .^ (-7+pH_vals);
% keq.gapdh.equilibrator
eQvals =    [1.7E-3;
            2.3E-3;
            4.0E-3;
            8.0E-3;
            1.62E-2;
            3.46E-2;
            6.56E-2;
            1.16E-1;
            1.78E-1;
            2.45E-1];
keq_gapdh_eQ = 0.05 * (eQvals);

% keq.pgk.eQuilibrator (forward reaction). Obtained in eQuilibrator
keq_pgk = [1/(7.4E-4); 1/(7.2E-4); 1/(6.9E-4); 1/(6.5E-4); 1/(6.1E-4); 1/(5.7E-4); 1/(5.5E-4); 1/(5.3E-4); 1/(5.2E-4); 1/(5.2E-4)];

figure(1)
% keq_gapdh
subplot(1,2,1)
plot(pH_vals,keq_gapdh_eQ,'.-','LineWidth',1.5,'MarkerSize',15)
hold on
plot(pH_vals,keq_gapdh,'.-','LineWidth',1.5,'MarkerSize',15)
legend('eQuilibrator','calculated','Location','northwest')
title('k_{eq,GAPDH} vs pH')
% keq_pgk
subplot(1,2,2)
plot(pH_vals,keq_pgk,'.-','LineWidth',1.5,'MarkerSize',15)
legend('eQuilibrator','Location','northwest')
title('k_{eq,PGK} vs pH')
