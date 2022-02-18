% % temp_keq_gapdh_check
% Reproducing the code that Bas sent on 12/05 to understand the difference
% with my simulations

% initial setup and parameters
% p = struct;
% % % % p.PGK_Vm = 1300/1;%60;                    % 21.7742
% % % % p.PGK_Keq = 1/(5.7E-4); % @7.06, but not much change %1.3514e+03
% % % % p.TDH1_Vmr = 0.0015*1000;                % 5.0631e-4
% % % % p.TDH1_Knad = 2.92;                 % 2.8799
% % % % p.TDH1_Knadh = 0.022;               % 0.0149
% % % % p.TDH1_Kbpg = 1.18;                 % 0.0437
% % % % p.TDH1_Kgap = 2.48;                 % 2.323
% % % % p.GAPDH_Keq = 3.46E-2 * 0.05;       % 8.5000e-05
p.PGK_Vm = 21.7742;
p.PGK_Keq = 1/(5.7E-4); % @7.06, but not much change %1.3514e+03
p.TDH1_Vmr = 0.0015;
p.TDH1_Knad = 2.8799;
p.TDH1_Knadh = 0.0149;
p.TDH1_Kbpg = 0.0437;
p.TDH1_Kgap = 2.323;
p.GAPDH_Keq = 0.0017;       % 8.5000e-05
% % Keq_vals = 0.05 * [1.7E-3, 2.3E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1, 1.78E-1, 2.45E-1];
% 50
% % % % % % % % Keq_vals = 50 * [1.7E-3, 2.3E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1, 1.78E-1, 2.45E-1];
    % Laura's gapdhrev vmax values
    %     0.0017
    %     0.0017
    %     0.0018
    %     0.0020
    %     0.0016
    %     0.0015 @pH 7.06
    %     0.0012
    %     0.0011
    %     0.0011
    %     0.0010
y0 = zeros(1,8);
% % % %     y0(1) = 5e-3; % P3G = y0(1); % start at 5e-3 M
% % % %     y0(2) = 1e-3; % ATP = y0(2); % start at 1e-3 M
% % % %     y0(3) = 0; % BPG = y0(3); % drops out
% % % %     y0(4) = 0; % ADP = y0(4); % start at 0
% % % %     y0(5) = 0; % NAD = y0(5); % start at 0
% % % %     y0(6) = 0; % GAP = y0(6); % start at 0
% % % %     y0(7) = 50e-3; % PHOS = y0(7); % start at 50e-3 M
% % % %     y0(8) = 0.15e-3; % NADH = y0(8); % start at 0.15e-3 M  
    y0(1) = 5; % P3G = y0(1); % start at 5e-3 M
    y0(2) = 1; % ATP = y0(2); % start at 1e-3 M
    y0(3) = 0; % BPG = y0(3); % drops out
    y0(4) = 0; % ADP = y0(4); % start at 0
    y0(5) = 0; % NAD = y0(5); % start at 0
    y0(6) = 0; % GAP = y0(6); % start at 0
    y0(7) = 500; % PHOS = y0(7); % start at 50e-3 M
    y0(8) = 0.0865; % NADH = y0(8); % start at 0.15e-3 M  
        
% % % % tspan = [0 1000];
tspan = [0 300];
options = [];

% % % % % % % % % simulation
% % % % % % % % n = length(Keq_vals);
% % % % % % % % simResults = cell(1,n); 
% % % % % % % % for j = 1:n
% % % % % % % %     p.GAPDH_Keq = Keq_vals(j);
% % % % % % % %     [t,y] = ode15s(@ode_gapdhr,tspan,y0,options,p);
% % % % % % % %     simResult{j}.t = t;
% % % % % % % %     simResult{j}.y = y;
% % % % % % % % end
% % % % % % % % 
% % % % % % % % % plot results
% % % % % % % % legendNames = {'p3g','atp','bpg','adp','nad','gap','phos','nadh'};
% % % % % % % % % figure(10001)
% % % % % % % % % for j = 1:n
% % % % % % % % %     for i = 1:length(y0)
% % % % % % % % % %     for i = 8
% % % % % % % % %         subplot(3,3,i)
% % % % % % % % %         plot(simResult{j}.t, simResult{j}.y(:,i),'.-')
% % % % % % % % %         hold on
% % % % % % % % %         legend(legendNames{i})
% % % % % % % % %     end
% % % % % % % % % end
% % % % % % % % figure(10002)
% % % % % % % % for j = 1:n
% % % % % % % %     for i = 8
% % % % % % % %         plot(simResult{j}.t, simResult{j}.y(:,i),'.-')
% % % % % % % %         hold on
% % % % % % % %         legend(legendNames{i})
% % % % % % % %     end
% % % % % % % % end


    
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
%     % values from Bas's model
%     
%     % parameters (pgk_keq is reversed between mine and Bas' simulation)
    p.PGK_Vm = 100;%DIF%21.7742;
    p.PGK_Keq = 0.00058;%DIF%1/(5.7E-4); % @7.06, but not much change %1.3514e+03
%     % DIF kinetics are not MM this time
    p.TDH1_Vmr = 0.03;%0.0015;
    p.TDH1_Knad = 2.8799;%
    p.TDH1_Knadh = 0.0149;%
    p.TDH1_Kbpg = 1;%0.0437;
    p.TDH1_Kgap = 2.323;%
%     p.GAPDH_Keq = 3.03836e-5;%0.0017;       % 8.5000e-05
    % initial concentrations
    y0 = zeros(1,8);
    y0(1) = 5e0; % P3G
    y0(2) = 1e0; % ATP
    y0(3) = 1e-3; % BPG
    y0(4) = 1e-5; % ADP
    y0(5) = 1e-5; % NAD
    y0(6) = 1e-5; % GAP
    y0(7) = 0; % PHOS
    y0(8) = 1.5e-1; % NADH
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    [t,y] = ode15s(@ode_gapdhr,tspan,y0,options,p);
    figure(1000), plot(t, y(:,8),'.-')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % FUNCTIONS CALLED: % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function v = ode_gapdhr(tspan,y0,p)

% recall
P3G = y0(1); % start at 5e-3 M      5
ATP = y0(2); % start at 1e-3 M      1
BPG = y0(3); % drops out            0
ADP = y0(4); % start at 0           0
NAD = y0(5); % start at 0           0
GAP = y0(6); % start at 0           0
PHOS = y0(7); % start at 50e-3 M    500
NADH = y0(8); % start at 0.15e-3 M  0.0865

% rates
v_PGK = p.PGK_Vm .* (BPG .* ADP - P3G .* ATP ./ p.PGK_Keq);
v_GAPDHrev = (p.TDH1_Vmr .* (BPG .* NADH - GAP .* NAD .* p.GAPDH_Keq)./(p.TDH1_Kbpg .* p.TDH1_Knadh))./...
    ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./...
    p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));

% odes
v(1) = + v_PGK; %P3G
v(2) = + v_PGK; %ATP
v(3) = - v_PGK - v_GAPDHrev; %BPG
v(4) = - v_PGK; %ADP
v(5) = + v_GAPDHrev; %NAD
v(6) = + v_GAPDHrev; %GAP
v(7) = + v_GAPDHrev; %PHOS
v(8) = - v_GAPDHrev; % NADH
v=v';

end


