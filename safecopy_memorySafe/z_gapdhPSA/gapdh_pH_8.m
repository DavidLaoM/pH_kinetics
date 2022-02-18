% % gapdh_ph_8
% start from toolbox implementation and go backwards.
% effect in change: 
    % gapdh.vf

% initial set of concentrations
y0 = zeros(1,7);
y0(1) = 5e0; % P3G = y0(1);
y0(2) = 1e0; % ATP = y0(2);
y0(3) = 0;%1e-3; % BPG = y0(3);
y0(4) = 0;%1e-5; % ADP = y0(4);
y0(5) = 0;%1e-5; % NAD = y0(5);
y0(6) = 0;%1e-5; % GAP = y0(6);
% y0(7) = 1.5e-1; % NADH = y0(7);
y0(7) = 9.84e-2;

% parameters
% p.gapdh.Vf = 0.03; % set as no enzyme gapdhrev
% % % % % p.gapdh.Vf = 0.0015; % % % %
% % p.gapdh.Kgap = 2.48;
% p.gapdh.Kgap = 2.3230;
% % p.gapdh.Knad = 2.92;
% p.gapdh.Knad = 2.8799;
% % p.gapdh.Keq_fwd = 0.0516 * 0.05;% p.gapdh.Keq_fwd = 3.04e-5;
% % % % % p.gapdh.Keq_fwd = 0.0017;
% p.gapdh.Keq_rev = 1/(0.0516 * 0.05);%3.04e-5;
% % p.gapdh.Kbpg = 1;
% p.gapdh.Kbpg = 0.0437;
% % p.gapdh.Knadh = 0.022;
% p.gapdh.Knadh = 0.0149;
% % p.pgk.Vf = 100;
% p.pgk.Vf = 21.7742;
% p.pgk.Kp3g = 0.53;
% p.pgk.Katp = 0.03;
% p.pgk.Keq_fwd = 1.7E3;
% p.pgk.Keq_rev = 5.8E-4;
% p.pgk.Kbpg = 0.003;
% p.pgk.Kadp = 0.2;
p.gapdh.Vf = 0.0033;
p.gapdh.Keq_fwd = 0.0017;
p.gapdh.Kgap = 2.3230;
p.gapdh.Knad = 2.8799;
p.gapdh.Kbpg = 0.0437;
p.gapdh.Knadh = 0.0149;
p.pgk.Vf = 21.7742;
p.pgk.Keq_fwd = 1.7E3;
%     p.gapdh.Vf = 0.03;

p.gapdh.Keq_fwd_array = 0.05 * [1.7E-3, 2.3E-3, 2.9E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1, 1.78E-1, 2.45E-1, 3.05E-1];
p.pgk.Keq_fwd_array = [1/(7.4E-4),    1/(7.2E-4),    1/(7.1E-4),    1/(6.9E-4),    1/(6.5E-4),    1/(6.1E-4),    1/(5.7E-4),    1/(5.5E-4),    1/(5.3E-4),    1/(5.2E-4),   1/(5.2E-4), 1/(5.1E-4)];
% p.gapdh.Keq_fwd_array = p.gapdh.Keq_fwd;
% p.pgk.Keq_fwd_array = p.pgk.Keq_fwd;

% simulation + visualization
tspan = [0 300];
options = [];
m = length(p.gapdh.Keq_fwd_array);
simResult = cell(1,m);
for j = 1:m
    p.gapdh.Keq_fwd = p.gapdh.Keq_fwd_array(j);
    p.pgk.Keq_fwd = p.pgk.Keq_fwd_array(j);
    [t,y] = ode15s(@ode_gapdh_pH_8,tspan,y0,options,p);
    r = calcRates_gapdh_pH_8(y,p);
    simResult{j}.t = t;
    simResult{j}.y = y;
    simResult{j}.r = r;
end
metNames = {'p3g','atp','bpg','adp','nad','gap','nadh'};
ratesNames = {'v_{pgk}','v_{gapdh.rev}'};
metIdxs = [7 9 4 6 2 1 5];
ratesIdxs = [4 1];

n = length(y0);
figure(101)
% figure(301)
for i = 1:n
    subplot(3,3,metIdxs(i))
    for j = 1:m
        if j == 6
            plot(simResult{j}.t,simResult{j}.y(:,i),'k.-')
        else
            plot(simResult{j}.t,simResult{j}.y(:,i),'-','color',[.5 .5 .5])
        end
        hold on
    end
    title(metNames{i})
end
suptitle('concentrations')
% set(101,'color','white')

figure(201)
% figure(401)

subplot(3,3,ratesIdxs(1)), 
for j = 1:m
    if j == 6
        plot(simResult{j}.t,simResult{j}.r(:,1),'k.-')
    else
        plot(simResult{j}.t,simResult{j}.r(:,1),'-','color',[.5 .5 .5])
    end
    hold on
end
title(ratesNames{1})

subplot(3,3,ratesIdxs(2)), 
for j = 1:m
    if j == 6
        plot(simResult{j}.t,simResult{j}.r(:,2),'k.-')
    else
        plot(simResult{j}.t,simResult{j}.r(:,2),'-','color',[.5 .5 .5])
    end
    hold on
end
title(ratesNames{2})

suptitle('reaction rates')
% set(201,'color','white')

% model
function v = ode_gapdh_pH_8(tspan,y0,p)

% recall initial concentrations
P3G = y0(1);
ATP = y0(2);
BPG = y0(3);
ADP = y0(4);
NAD = y0(5);
GAP = y0(6);
NADH = y0(7);

% rates
% % Topology manuscript
vPGK = p.pgk.Vf .* (BPG .* ADP - P3G .* ATP ./ p.pgk.Keq_fwd);
vGAPDHrev = (p.gapdh.Vf .* ...
    (BPG .* NADH - GAP .* NAD ./ (1./ p.gapdh.Keq_fwd)) .* ...
    (1./(p.gapdh.Kbpg .* p.gapdh.Knadh)))./...
    ((1 + NAD ./ p.gapdh.Knad + NADH ./ p.gapdh.Knadh) .* (1 + BPG ./ p.gapdh.Kbpg + GAP ./ p.gapdh.Kgap));
% % vGAPDHrev = -(p.gapdh.Vf .* ...
% %     (BPG .* NADH .* p.gapdh.Keq_fwd - GAP .* NAD ./ (1./ 1)) .* ...
% %     (1./(p.gapdh.Kbpg .* p.gapdh.Knadh)))./...
% %     ((1 + NAD ./ p.gapdh.Knad + NADH ./ p.gapdh.Knadh) .* (1 + BPG ./ p.gapdh.Kbpg + GAP ./ p.gapdh.Kgap));
%     ((1 + GAP ./ p.gapdh.Kgap + NAD ./ p.gapdh.Knad) .* (1 + BPG ./ p.gapdh.Kbpg + NADH ./ p.gapdh.Knadh));
% vGAPDHrev = -p.gapdh.Vf .* ...
%     (1./(p.gapdh.Kgap .* p.gapdh.Knad))
%     (GAP .* NAD - BPG .* NADH ./ p.gapdh.Keq_fwd) ./ ...
%     ((1 + GAP ./ p.gapdh.Kgap + NAD ./ p.gapdh.Knad) .* (1 + BPG ./ p.gapdh.Kbpg + NADH ./ p.gapdh.Knadh));
% %     ((1 + NAD ./ p.gapdh.Knad + NADH ./ p.gapdh.Knadh) .* (1 + BPG ./ p.gapdh.Kbpg + GAP ./ p.gapdh.Kgap));

% % % % % % Topology under study
% % % % vPGKrev = p.pgk.Vf .* ...
% % % %     (1 ./ (p.pgk.Kp3g .* p.pgk.Katp)) .* ...
% % % %     (P3G .* ATP - BPG .* ADP .* p.pgk.Keq_fwd) ./ ...
% % % %     ((1 + P3G ./ p.pgk.Kp3g + ATP ./ p.pgk.Katp) .* (1 + BPG ./ p.pgk.Kbpg + ADP ./ p.pgk.Kadp));
% % % % vPGK = - vPGKrev;
% vGAPDH = p.gapdh.Vf .* ... &1.0274e+04 .* ... %(70/1184) .* ...
%     (GAP .* NAD .* p.gapdh.Keq_fwd - BPG .* NADH ./ 1) .* ...
%     (1 ./ (p.gapdh.Kgap .* p.gapdh.Knad)) ./ ...
%     ((1 + GAP ./ p.gapdh.Kgap + NAD ./ p.gapdh.Knad) .* (1 + BPG ./ p.gapdh.Kbpg + NADH ./ p.gapdh.Knadh));
% %     ((1 + NAD ./ p.gapdh.Knad + NADH ./ p.gapdh.Knadh) .* (1 + BPG ./ p.gapdh.Kbpg + GAP ./ p.gapdh.Kgap));
% vGAPDHrev = -vGAPDH;

vGAPDH = p.gapdh.Vf .* ... &1.0274e+04 .* ... %(70/1184) .* ...
    (GAP .* NAD - BPG .* NADH ./ p.gapdh.Keq_fwd) .* ...
    (1 ./ (p.gapdh.Kgap .* p.gapdh.Knad)) ./ ...
    ((1 + GAP ./ p.gapdh.Kgap + NAD ./ p.gapdh.Knad) .* (1 + BPG ./ p.gapdh.Kbpg + NADH ./ p.gapdh.Knadh));
%     ((1 + NAD ./ p.gapdh.Knad + NADH ./ p.gapdh.Knadh) .* (1 + BPG ./ p.gapdh.Kbpg + GAP ./ p.gapdh.Kgap));
vGAPDHrev = -vGAPDH;

% mass balance
v(1) = + vPGK; % p3g
v(2) = + vPGK; % atp
v(3) = - vPGK - vGAPDHrev; % bpg
v(4) = - vPGK; % adp
v(5) = + vGAPDHrev; % nad
v(6) = + vGAPDHrev; % gap
v(7) = - vGAPDHrev; % nadh
% v(1) = - vPGKrev; % p3g
% v(2) = - vPGKrev; % atp
% v(3) = + vPGKrev + vGAPDH; % bpg
% v(4) = + vPGKrev; % adp
% v(5) = - vGAPDH; % nad
% v(6) = - vGAPDH; % gap
% v(7) = + vGAPDH; % nadh
v = v';    

end

% calcRates
function r = calcRates_gapdh_pH_8(y,p)

% recall initial concentrations
P3G = y(:,1);
ATP = y(:,2);
BPG = y(:,3);
ADP = y(:,4);
NAD = y(:,5);
GAP = y(:,6);
NADH = y(:,7);

% calcRates
% vGAPDH = p.gapdh.Vf .*...
%     (1./(p.gapdh.Kgap .* p.gapdh.Knad)) .* ...
%     (GAP .* NAD - BPG .* NADH ./ p.gapdh.Keq_fwd) .* ...
%     ((1 + GAP./p.gapdh.Kgap + NAD./p.gapdh.Knad) .* (1 + BPG./p.gapdh.Kbpg + NADH./p.gapdh.Knadh));
% % vGAPDHrev = p.gapdh.Vf .*...
% %     (1./(p.gapdh.Kbpg .* p.gapdh.Knadh)) .* ...
% %     (BPG .* NADH - GAP .* NAD .* p.gapdh.Keq_fwd) .* ...
% %     ((1 + NADH./p.gapdh.Knadh + NAD./p.gapdh.Knad) .* (1 + BPG./p.gapdh.Kbpg + GAP./p.gapdh.Kgap));
% vPGK = p.pgk.Vf .*...
%     (1./(p.pgk.Kbpg .* p.pgk.Kadp)) .* ...
%     (BPG .* ADP - P3G .* ATP ./ p.pgk.Keq_fwd) .* ...
%     ((1 + P3G./p.pgk.Kp3g + ATP./p.pgk.Katp) .* (1 + BPG./p.pgk.Kbpg + ADP./p.pgk.Kadp));
% % vPGK = p.pgk.Vf .* (BPG .* ADP - P3G .* ATP ./ p.pgk.Keq_fwd);
% r = [vPGK, vGAPDH];
vPGK = p.pgk.Vf .* (BPG .* ADP - P3G .* ATP ./ p.pgk.Keq_fwd);
vGAPDHrev = (p.gapdh.Vf .* (BPG .* NADH - GAP .* NAD .* p.gapdh.Keq_fwd)./(p.gapdh.Kbpg .* p.gapdh.Knadh))./...
    ((1 + NAD ./ p.gapdh.Knad + NADH ./ p.gapdh.Knadh) .* (1 + BPG ./...
    p.gapdh.Kbpg + GAP ./ p.gapdh.Kgap));
r = [vPGK, -vGAPDHrev];

end

