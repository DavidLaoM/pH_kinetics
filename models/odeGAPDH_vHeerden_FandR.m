function [v] = odeGAPDH_vHeerden_FandR(tspan,y0,p,f,data,setup)
%ODEGAPDHR Summary of this function goes here
%   Kinetics from vanHeerden 2014
%   Order of events
%       Select initial concentraitons
%       Calculate v
%       Mass balances
ode_pH = setup.ode_pH;
typeVm = setup.typeVm;

% select initial points
% reverse
rP3G = y0(1);
rATP = y0(2);
rBPG = y0(3);
rADP = y0(4);
rNAD = y0(5);
rGAP = y0(6);
rPHOS = y0(7);
rNADH = y0(8);
% forward
fP3G = y0(9);
fATP = y0(10);
fBPG = y0(11);
fADP = y0(12);
fNAD = y0(13);
fGAP = y0(14);
fPHOS = y0(15);
fNADH = y0(16);

% calulate v (rateEquations)

switch typeVm
    case 'specific'
        switch ode_pH
            case 'on'
                H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
                v_GAPDHr = (-(p.TDH1_Vmr_R .* rBPG .* rNADH .* H ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf_R .* rNAD .* rGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap));
                v_GAPDHf = (-(p.TDH1_Vmr_F .* fBPG .* fNADH .* H ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf_F .* fNAD .* fGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + fNAD ./ p.TDH1_Knad + fNADH ./ p.TDH1_Knadh) .* (1 + fBPG ./ p.TDH1_Kbpg + fGAP ./ p.TDH1_Kgap));
            case 'on_vmf_vmr'
                H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
                v_GAPDHr = -p.TDH1_Vmr_R .* rBPG .* rNADH .* H + p.TDH1_Vmf_R .* rNAD .* rGAP;
                v_GAPDHf = -p.TDH1_Vmr_F .* fBPG .* fNADH .* H + p.TDH1_Vmf_F .* fNAD .* fGAP;
            case 'on_vm_keq' %now 'p.TDH1_Vmr_R' is a dummy for Keq. Boundaries for x(6) may have to be increased
                H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
                v_GAPDHr = p.TDH1_Vmf_R .* ( rNAD .* rGAP - p.TDH1_Vmr_R .* rBPG .* rNADH .* H);
                v_GAPDHf = p.TDH1_Vmf_F .* ( rNAD .* rGAP - p.TDH1_Vmr_F .* rBPG .* rNADH .* H);
            case 'on_revMM_bitri'      
                H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
                v_GAPDHr = p.TDH1_Vmf_R .* ( rNAD .* rGAP - p.TDH1_Vmr_R .* rBPG .* rNADH .* H)./((p.TDH1_Kgap .* p.TDH1_Knad).*((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap)));
                v_GAPDHf = p.TDH1_Vmf_F .* ( rNAD .* rGAP - p.TDH1_Keq .* rBPG .* rNADH .* H)./((p.TDH1_Kgap .* p.TDH1_Knad).*((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap)));
            otherwise
                v_GAPDHr = (-(p.TDH1_Vmr_R .* rBPG .* rNADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf_R .* rNAD .* rGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap));
                v_GAPDHf = (-(p.TDH1_Vmr_F .* fBPG .* fNADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf_F .* fNAD .* fGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + fNAD ./ p.TDH1_Knad + fNADH ./ p.TDH1_Knadh) .* (1 + fBPG ./ p.TDH1_Kbpg + fGAP ./ p.TDH1_Kgap));
        end   
    case 'common'
        switch ode_pH
            case 'on'
                H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
                v_GAPDHr = (-(p.TDH1_Vmr .* rBPG .* rNADH .* H ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* rNAD .* rGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap));
                v_GAPDHf = (-(p.TDH1_Vmr .* fBPG .* fNADH .* H ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* fNAD .* fGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + fNAD ./ p.TDH1_Knad + fNADH ./ p.TDH1_Knadh) .* (1 + fBPG ./ p.TDH1_Kbpg + fGAP ./ p.TDH1_Kgap));
            case 'on_vmf_vmr'
                H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
                v_GAPDHr = -p.TDH1_Vmr .* rBPG .* rNADH .* H + p.TDH1_Vmf .* rNAD .* rGAP;
                v_GAPDHf = -p.TDH1_Vmr .* fBPG .* fNADH .* H + p.TDH1_Vmf .* fNAD .* fGAP;
            case 'on_vm_keq' %now 'p.TDH1_Vmr_R' is a dummy for Keq. Boundaries for x(6) may have to be increased
                H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
                v_GAPDHr = p.TDH1_Vmf .* ( rNAD .* rGAP - p.TDH1_Vmr .* rBPG .* rNADH .* H);
                v_GAPDHf = p.TDH1_Vmf .* ( rNAD .* rGAP - p.TDH1_Vmr .* rBPG .* rNADH .* H);
            case 'on_revMM_bitri'
                H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
                v_GAPDHr = p.TDH1_Vmf .* ( rNAD .* rGAP - p.TDH1_Vmr .* rBPG .* rNADH .* H)./((p.TDH1_Kgap .* p.TDH1_Knad).*((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap)));
                v_GAPDHf = p.TDH1_Vmf .* ( rNAD .* rGAP - p.TDH1_Vmr .* rBPG .* rNADH .* H)./((p.TDH1_Kgap .* p.TDH1_Knad).*((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap)));
            otherwise
                v_GAPDHr = (-(p.TDH1_Vmr .* rBPG .* rNADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* rNAD .* rGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap));
                v_GAPDHf = (-(p.TDH1_Vmr .* fBPG .* fNADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* fNAD .* fGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + fNAD ./ p.TDH1_Knad + fNADH ./ p.TDH1_Knadh) .* (1 + fBPG ./ p.TDH1_Kbpg + fGAP ./ p.TDH1_Kgap));
        end      
    otherwise
        disp('No typeVm is selected.');
end  
v_PGKr = p.PGK_Vm .* (rBPG .* rADP - rP3G .* rATP ./ p.PGK_Keq);
v_PGKf = p.PGK_Vm .* (fBPG .* fADP - fP3G .* fATP ./ p.PGK_Keq);
% v_PGK = p.PGK_Vm .* (BPG .* ADP - P3G .* ATP ./ p.PGK_Keq) ./ (p.PGK_Katp .* p.PGK_Kp3g .* (1 + BPG ./ p.PGK_Kbpg + P3G ./ p.PGK_Kp3g) .* (1 + ADP ./ p.PGK_Kadp + ATP ./ p.PGK_Katp));

% mass balances (assumed ideally mixed batch system configuration)
% reverse
v(1) = + v_PGKr; %P3G
v(2) = + v_PGKr; %ATP
v(3) = - v_PGKr + v_GAPDHr; %BPG
v(4) = - v_PGKr; %ADP
v(5) = - v_GAPDHr; %NAD
v(6) = - v_GAPDHr; %GAP
v(7) = - v_GAPDHr; %PHOS
v(8) = + v_GAPDHr; % NADH
% forward
v(9) = + v_PGKf; %P3G
v(10) = + v_PGKf; %ATP
v(11) = - v_PGKf + v_GAPDHf; %BPG
v(12) = - v_PGKf; %ADP
v(13) = - v_GAPDHf; %NAD
v(14) = - v_GAPDHf; %GAP
v(15) = - v_GAPDHf; %PHOS
v(16) = + v_GAPDHf; % NADH
v=v';

end
%
