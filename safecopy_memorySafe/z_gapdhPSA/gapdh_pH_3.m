% % gapdh_ph_3
% pgk_rev

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
p.gapdhrev.Vf = 0;%0.03; % set as no enzyme gapdhrev
p.gapdhrev.Ks1 = 2.48;
p.gapdhrev.Ks2 = 2.92;
p.gapdhrev.Keq = 3.04e-5;
p.gapdhrev.pH = 6.8;
p.gapdhrev.Kp1 = 1;
p.gapdhrev.Kp2 = 0.022;

p.pgk.Vf = 100;
p.pgk.Kp3g = 0.53;
p.pgk.Katp = 0.03;
p.pgk.Keq_fwd = 1.7E3;
p.pgk.Keq_rev = 5.8E-4;
p.pgk.Kbpg = 0.003;
p.pgk.Kadp = 0.2;

% simulation + visualization
tspan = [0 300];
options = [];
[t,y] = ode15s(@ode_gapdh_pH_3,tspan,y0,options,p);
r = calcRates_gapdh_pH_3(y,p);
metNames = {'p3g','atp','bpg','adp','nad','gap','nadh'};
ratesNames = {'v_{pgk}','v_{gapdh.rev}'};
metIdxs = [7 9 4 6 2 1 5];
ratesIdxs = [4 1];

n = length(y0);
figure(101)
% figure(301)
for i = 1:n
    subplot(3,3,metIdxs(i))
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
subplot(3,3,ratesIdxs(1)), plot(t,r(:,1),'k.-'), title(ratesNames{1})
subplot(3,3,ratesIdxs(2)), plot(t,r(:,2),'k.-'), title(ratesNames{2})
suptitle('reaction rates')
% set(201,'color','white')

% model
function v = ode_gapdh_pH_3(tspan,y0,p)

% recall initial concentrations
P3G = y0(1);
ATP = y0(2);
BPG = y0(3);
ADP = y0(4);
NAD = y0(5);
GAP = y0(6);
NADH = y0(7);

% rates
vGAPDH = 0;
% vPGK = p.pgk.Vf .*...
%     (1./(p.pgk.Kbpg .* p.pgk.Kadp)) .* ...
%     (BPG .* ADP - P3G .* ATP ./ p.pgk.Keq_fwd) .* ...
%     ((1 + P3G./p.pgk.Kp3g + ATP./p.pgk.Katp) .* (1 + BPG./p.pgk.Kbpg + ADP./p.pgk.Kadp));
vPGKrev = p.pgk.Vf .*...
    (1./(p.pgk.Kp3g .* p.pgk.Katp)) .* ...
    (P3G .* ATP - BPG .* ADP./ p.pgk.Keq_rev) ./ ...
    ((1 + P3G./p.pgk.Kp3g + ATP./p.pgk.Katp) .* (1 + BPG./p.pgk.Kbpg + ADP./p.pgk.Kadp));

% mass balance
v(1) = - vPGKrev; % p3g
v(2) = - vPGKrev; % atp
v(3) = + vPGKrev + vGAPDH; % bpg
v(4) = + vPGKrev; % adp
v(5) = - vGAPDH; % nad
v(6) = - vGAPDH; % gap
v(7) = + vGAPDH; % nadh
v = v';    

end

% calcRates
function r = calcRates_gapdh_pH_3(y,p)
% recall initial concentrations
P3G = y(:,1);
ATP = y(:,2);
BPG = y(:,3);
ADP = y(:,4);
NAD = y(:,5);
GAP = y(:,6);
NADH = y(:,7);
% calcRates
vGAPDH = zeros(size(y(:,1)));
% vPGK = p.pgk.Vf .*...
%     (1./(p.pgk.Kbpg .* p.pgk.Kadp)) .* ...
%     (BPG .* ADP - P3G .* ATP ./ p.pgk.Keq_fwd) .* ...
%     ((1 + P3G./p.pgk.Kp3g + ATP./p.pgk.Katp) .* (1 + BPG./p.pgk.Kbpg + ADP./p.pgk.Kadp));
vPGKrev = p.pgk.Vf .*...
    (1./(p.pgk.Kp3g .* p.pgk.Katp)) .* ...
    (P3G .* ATP - BPG .* ADP./ p.pgk.Keq_rev) ./ ...
    ((1 + P3G./p.pgk.Kp3g + ATP./p.pgk.Katp) .* (1 + BPG./p.pgk.Kbpg + ADP./p.pgk.Kadp));
r = [-vPGKrev, vGAPDH];
end

