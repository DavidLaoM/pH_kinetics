function [simResult] = simSys_gapdhFWDREV(xtemp,data,setup)
% file where simulation is run
enzymeName = setup.enzymeName;
constantVm = setup.constantVm;
ode = setup.ode;
plotEachSim = setup.plotEachSim;
legendamets = setup.PSAmets;
sourceVm = setup.sourceVm;
% disp(data.i); 

% set params and ode
odefun = @odeGAPDH_fwd_rev;
p.TDH1_Kgap=10.^xtemp(1).*2.48; % mM
p.TDH1_Kbpg=10.^xtemp(2).*1.18; % mM
p.TDH1_Knad=10.^xtemp(3).*2.92; %mM
p.TDH1_Knadh=10.^xtemp(4).*0.022; % mM
p.TDH1_Vmf = data.chosenVmaxFWD * 10 .^ xtemp(5); %(UNIT!) 
p.TDH1_Vmr = data.chosenVmaxREV * 10 .^ xtemp(6); %(UNIT!) 
p.GAPDH_Keq = data.chosenKeqGAPDH;
p.PGK_Vm = 1306.45 / 60;
p.PGK_Keq = data.chosenKeqPGK;

% initial concentrations
if setup.directionFWD == 1
    P3Go = 0;
    ATPo = 0;
    BPGo = 0;
    ADPo = 10;
    NADo = 1;
    GAPo = 5.8;
    PHOSo = 0;
    NADHo = data.chosenNADini;
elseif setup.directionREV == 1
    P3Go = 5;
    ATPo = 1;
    BPGo = 0;
    ADPo = 0;
    NADo = 0;
    GAPo = 0;
    PHOSo = 500;
    NADHo = data.chosenNADini;
else
    disp('Warning: no direction (fwd or rev) has been selected.')
end

% reaction setup and run
tspan   = [0 600];  % time [300 s (6 min)]
y0 = [P3Go ATPo BPGo ADPo NADo GAPo PHOSo NADHo];
options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
% run ode15s
f = 1; % not being used for PSA
[t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);
simResult.t = t;
simResult.y = y;

end


