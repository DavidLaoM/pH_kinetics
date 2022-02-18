% % gapdh_ph_1
% topology Bas

% initial set of concentrations
y0 = zeros(1,7);
y0(1) = 5e0; % P3G = y0(1);
y0(2) = 1e0; % ATP = y0(2);
y0(3) = 1e-3; % BPG = y0(3);
y0(4) = 1e-5; % ADP = y0(4);
y0(5) = 1e-5; % NAD = y0(5);
y0(6) = 1e-5; % GAP = y0(6);
y0(7) = 1.5e-1; % NADH = y0(7);

% parameters
p.gapdhrev.Vf = 0.03;
p.gapdhrev.Ks1 = 2.48;
p.gapdhrev.Ks2 = 2.92;
p.gapdhrev.Keq = 3.04e-5;
p.gapdhrev.pH = 6.8;
p.gapdhrev.Kp1 = 1;
p.gapdhrev.Kp2 = 0.022;

p.pgk.Vf = 100;
p.pgk.Ks1 = 0.53;
p.pgk.Ks2 = 0.03;
p.pgk.Keq = 0.00058;
p.pgk.Kp1 = 0.003;
p.pgk.Kp2 = 0.2;

% simulation + visualization
tspan = [0 300];
options = [];
[t,y] = ode15s(@ode_gapdh_pH_1,tspan,y0,options,p);
r = calcRates_gapdh_pH_1(y,p);
metNames = {'p3g','atp','bpg','adp','nad','gap','nadh'};
ratesNames = {'v_{pgk}','v_{gapdh.rev}'};

n = length(y0);
figure(101)
% figure(301)
for i = 1:n
    subplot(3,3,i)
    if i == n
        plot(t,y(:,i),'k.-','LineWidth',2)
    else
        plot(t,y(:,i),'k.-')
    end
    title(metNames{i})
end
suptitle('concentrations')
% set(101,'color','white')

figure(201)
% figure(401)
subplot(3,3,1), plot(t,r(:,1),'k.-'), title(ratesNames{1})
subplot(3,3,2), plot(t,r(:,2),'k.-'), title(ratesNames{2})
suptitle('reaction rates')
% set(201,'color','white')

% model
function v = ode_gapdh_pH_1(tspan,y0,p)
% recall initial concentrations
P3G = y0(1);
ATP = y0(2);
BPG = y0(3);
ADP = y0(4);
NAD = y0(5);
GAP = y0(6);
NADH = y0(7);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%     % parameters (pgk_keq is reversed between mine and Bas' simulation)
p.PGK_Vm = 100;%DIF%21.7742;
% p.PGK_Keq = 1754;%0.00058;%DIF%1/(5.7E-4); % @7.06, but not much change %1.3514e+03
p.PGK_Keq = 0.00058;
%     % DIF kinetics are not MM this time
p.TDH1_Vmr = 0.03;%0.0015;
p.TDH1_Knad = 2.8799;%
p.TDH1_Knadh = 0.0149;%
p.TDH1_Kbpg = 1;%0.0437;
p.TDH1_Kgap = 2.323;%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% rates
v_PGK = p.PGK_Vm .* (BPG .* ADP - P3G .* ATP ./ p.PGK_Keq);
v_GAPDHrev = (p.TDH1_Vmr .* (BPG .* NADH - GAP .* NAD .* p.GAPDH_Keq)./(p.TDH1_Kbpg .* p.TDH1_Knadh))./...
    ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./...
    p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
vPGK = - v_PGK;
vGAPDHrev = - v_GAPDHrev;
% % % % % reaction rates simulation
% vPGK = p.PGK_Vm .* (P3G.*ATP - BPG.*ADP ./ p.PGK_Keq);
% % % % vPGK = p.pgk.Vf .* ...
% % % %     (1./(p.pgk.Ks1.*p.pgk.Ks2)) .* ...
% % % %     (P3G.*ATP - BPG.*ADP ./ p.pgk.Keq) .* ...
% % % %     ((1 + P3G./p.pgk.Ks1 + ATP./p.pgk.Ks2).*(1 + BPG./p.pgk.Kp1 + ADP./p.pgk.Kp2));
% % % % vGAPDHrev = p.gapdhrev.Vf .* ...
% % % %     (1./(p.gapdhrev.Ks1.*p.gapdhrev.Ks2)) .* ...
% % % %     (GAP.*NAD - BPG.*NADH ./ (p.gapdhrev.Keq + 10 .^ (-7 + p.gapdhrev.pH))) .* ...
% % % %     ((1 + GAP./p.gapdhrev.Ks1 + NAD./p.gapdhrev.Ks2) .* (1 + BPG./p.gapdhrev.Kp1 + NADH./p.gapdhrev.Kp2));
% vGAPDHrev = p.TDH1_Vmr .* ...
%     (1./(p.TDH1_Kbpg .* p.TDH1_Knadh)) .* ...
%     (GAP.*NAD - BPG.*NADH ./ (p.gapdhrev.Keq)) .* ...
%     ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));

% % odes
% v(1) = + v_PGK; %P3G
% v(2) = + v_PGK; %ATP
% v(3) = - v_PGK - v_GAPDHrev; %BPG
% v(4) = - v_PGK; %ADP
% v(5) = + v_GAPDHrev; %NAD
% v(6) = + v_GAPDHrev; %GAP
% v(7) = - v_GAPDHrev; % NADH
% v=v';
% mass balance
v(1) = - vPGK; % p3g
v(2) = - vPGK; % atp
v(3) = + vPGK + vGAPDHrev; % bpg
v(4) = + vPGK; % adp
v(5) = - vGAPDHrev; % nad
v(6) = - vGAPDHrev; % gap
v(7) = + vGAPDHrev; % nadh
v = v';    

end

% calcRates
function r = calcRates_gapdh_pH_1(y,p)
% recall initial concentrations
P3G = y(:,1);
ATP = y(:,2);
BPG = y(:,3);
ADP = y(:,4);
NAD = y(:,5);
GAP = y(:,6);
NADH = y(:,7);
% calcRates
vPGK = p.pgk.Vf .* ...
    (1./(p.pgk.Ks1.*p.pgk.Ks2)) .* ...
    (P3G.*ATP - BPG.*ADP ./ p.pgk.Keq) .* ...
    ((1 + P3G./p.pgk.Ks1 + ATP./p.pgk.Ks2).*(1 + BPG./p.pgk.Kp1 + ADP./p.pgk.Kp2));
vGAPDHrev = p.gapdhrev.Vf .* ...
    (1./(p.gapdhrev.Ks1.*p.gapdhrev.Ks2)) .* ...
    (GAP.*NAD - BPG.*NADH ./ (p.gapdhrev.Keq + 10 .^ (-7 + p.gapdhrev.pH))) .* ...
    ((1 + GAP./p.gapdhrev.Ks1 + NAD./p.gapdhrev.Ks2) .* (1 + BPG./p.gapdhrev.Kp1 + NADH./p.gapdhrev.Kp2));
r = [vPGK, vGAPDHrev];
end


% % rates
% v_PGK = p.pgk.Vf .* (BPG .* ADP - P3G .* ATP ./ p.pgk.Keq);
% v_GAPDHrev = (p.gapdhrev.Vf .* (BPG .* NADH - GAP .* NAD .* p.gapdhrev.Keq)./(p.gapdhrev.Ks1.*p.gapdhrev.Ks2))./...
%     ((1 + NAD ./ p.gapdhrev.Ks2 + NADH ./ p.gapdhrev.Kp2) .* (1 + BPG ./...
%     p.gapdhrev.Kp1 + GAP ./ p.gapdhrev.Ks1));
% 
% % odes
% v(1) = + v_PGK; %P3G
% v(2) = + v_PGK; %ATP
% v(3) = - v_PGK - v_GAPDHrev; %BPG
% v(4) = - v_PGK; %ADP
% v(5) = + v_GAPDHrev; %NAD
% v(6) = + v_GAPDHrev; %GAP
% v(7) = - v_GAPDHrev; % NADH
% v=v';

