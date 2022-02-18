function [vObs, vextra] = calcRates_FandR(xtemp,simResult,data,setup)
enzymeName = setup.enzymeName;
constantVm = setup.constantVm;
ode = setup.ode;
ode_pH = setup.ode_pH;
sourceVm = setup.sourceVm;
typeVm = setup.typeVm;

switch enzymeName
    case 'gapdhr'
        % metabolites
        % reverse
        rP3G = simResult.y(:,1);
        rATP = simResult.y(:,2);
        rBPG = simResult.y(:,3);
        rADP = simResult.y(:,4);
        rNAD = simResult.y(:,5);
        rGAP = simResult.y(:,6);
        rPHOS = simResult.y(:,7);
        rNADH = simResult.y(:,8);
        % forward
        fP3G = simResult.y(:,9);
        fATP = simResult.y(:,10);
        fBPG = simResult.y(:,11);
        fADP = simResult.y(:,12);
        fNAD = simResult.y(:,13);
        fGAP = simResult.y(:,14);
        fPHOS = simResult.y(:,15);
        fNADH = simResult.y(:,16);
        switch ode
            case 'vanHeerden2014'
                % GAPDH     
                p.TDH1_Kgap = 10.^xtemp(2).*2.48; % mM
                p.TDH1_Kbpg = 10.^xtemp(3).*1.18; % mM
                p.TDH1_Knad = 10.^xtemp(4).*2.92; %mM
                p.TDH1_Knadh = 10.^xtemp(5).*0.022; % mM
                p.TDH1_Keq = data.chosenKeqGAPDH; % []                
                switch typeVm
                    case 'specific'
                        switch sourceVm
                            case 'literature'
                                p.TDH1_Vmf_R = 10.^xtemp(1).*1184.52/60 / data.chosenDF;% mM s^{-1} 
                                p.TDH1_Vmr_R = 10.^xtemp(6).*6549.8/60 / data.chosenDF; % mM s^{-1}
                                p.TDH1_Vmf_F = 10.^xtemp(7).*1184.52/60 / data.chosenDF;% mM s^{-1} 
                                p.TDH1_Vmr_F = 10.^xtemp(8).*6549.8/60 / data.chosenDF; % mM s^{-1}
                            case 'experimentalSlopes'
                                p.TDH1_Vmf_R = 10.^xtemp(1).*setup.exp_vmax_gapdhf(data.i) / data.chosenDF;% mM s^{-1}
                                p.TDH1_Vmr_R = 10.^xtemp(6).*setup.exp_vmax_gapdhr(data.i) / data.chosenDF; % mM s^{-1} 
                                p.TDH1_Vmf_F = 10.^xtemp(7).*setup.exp_vmax_gapdhf(data.i) / data.chosenDF;% mM s^{-1}
                                p.TDH1_Vmr_F = 10.^xtemp(8).*setup.exp_vmax_gapdhr(data.i) / data.chosenDF; % mM s^{-1} 
                            case 'experimentalSlopesFixed'
                                p.TDH1_Vmf_R = 10.^xtemp(1).*setup.exp_vmax_gapdhf(6) / data.chosenDF;% mM s^{-1}
                                p.TDH1_Vmr_R = 10.^xtemp(6).*setup.exp_vmax_gapdhr(6) / data.chosenDF; % mM s^{-1}
                                p.TDH1_Vmf_F = 10.^xtemp(7).*setup.exp_vmax_gapdhf(6) / data.chosenDF;% mM s^{-1}
                                p.TDH1_Vmr_F = 10.^xtemp(8).*setup.exp_vmax_gapdhr(6) / data.chosenDF; % mM s^{-1}
                            otherwise
                                disp('No source for vmax has been selected');
                        end
                    case 'common'
                        switch sourceVm
                            case 'literature'
                                p.TDH1_Vmf = 10.^xtemp(1).*1184.52/60 / data.chosenDF;% mM s^{-1} 
                                p.TDH1_Vmr = 10.^xtemp(6).*6549.8/60 / data.chosenDF; % mM s^{-1}
                            case 'experimentalSlopes'
                                p.TDH1_Vmf = 10.^xtemp(1).*setup.exp_vmax_gapdhf(data.i) / data.chosenDF;% mM s^{-1}
                                p.TDH1_Vmr = 10.^xtemp(6).*setup.exp_vmax_gapdhr(data.i) / data.chosenDF; % mM s^{-1} %.*data.chosenVmax 
                            case 'experimentalSlopesFixed'
                                p.TDH1_Vmf = 10.^xtemp(1).*setup.exp_vmax_gapdhf(6) / data.chosenDF;% mM s^{-1}
                                p.TDH1_Vmr = 10.^xtemp(6).*setup.exp_vmax_gapdhr(6) / data.chosenDF; % mM s^{-1} %.*data.chosenVmax
                            otherwise
                                disp('No source for vmax has been selected');
                        end
                    otherwise
                        disp('No typeVm is selected.');
                end               
                % PGK
                p.PGK_Keq = data.chosenKeqPGK; % [] %/10
                p.PGK_Vm = 1306.45 / 60 * setup.excessPGK / data.chosenDF; % mM s^{-1} % corrected to make it appear in excess
                p.PGK_Katp = 0.3; % mM
                p.PGK_Kp3g = 0.53; % mM
                p.PGK_Kbpg = 0.003; % mM
                p.PGK_Kadp = 0.2; % mM
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
                                v_GAPDHf = p.TDH1_Vmf_F .* ( fNAD .* fGAP - p.TDH1_Vmr_F .* fBPG .* fNADH .* H);
                            case 'on_revMM_bitri'        
                                H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
                                v_GAPDHr = p.TDH1_Vmf_R .* ( rNAD .* rGAP - p.TDH1_Vmr_R .* rBPG .* rNADH .* H)./((p.TDH1_Kgap .* p.TDH1_Knad).*((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap)));
                                v_GAPDHf = p.TDH1_Vmf_F .* ( fNAD .* fGAP - p.TDH1_Vmr_F .* fBPG .* fNADH .* H)./((p.TDH1_Kgap .* p.TDH1_Knad).*((1 + fNAD ./ p.TDH1_Knad + fNADH ./ p.TDH1_Knadh) .* (1 + fBPG ./ p.TDH1_Kbpg + fGAP ./ p.TDH1_Kgap)));
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
                                v_GAPDHf = p.TDH1_Vmf .* ( fNAD .* fGAP - p.TDH1_Vmr .* fBPG .* fNADH .* H);
                            case 'on_revMM_bitri' 
                                H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
                                v_GAPDHr = p.TDH1_Vmf .* ( rNAD .* rGAP - p.TDH1_Vmr .* rBPG .* rNADH .* H)./((p.TDH1_Kgap .* p.TDH1_Knad).*((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap)));
                                v_GAPDHf = p.TDH1_Vmf .* ( fNAD .* fGAP - p.TDH1_Vmr .* fBPG .* fNADH .* H)./((p.TDH1_Kgap .* p.TDH1_Knad).*((1 + fNAD ./ p.TDH1_Knad + fNADH ./ p.TDH1_Knadh) .* (1 + fBPG ./ p.TDH1_Kbpg + fGAP ./ p.TDH1_Kgap)));
                            otherwise
                                v_GAPDHr = (-(p.TDH1_Vmr .* rBPG .* rNADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* rNAD .* rGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap));
                                v_GAPDHf = (-(p.TDH1_Vmr .* fBPG .* fNADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* fNAD .* fGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + fNAD ./ p.TDH1_Knad + fNADH ./ p.TDH1_Knadh) .* (1 + fBPG ./ p.TDH1_Kbpg + fGAP ./ p.TDH1_Kgap));
                        end      
                    otherwise
                        disp('No typeVm is selected.');
                end
                v_PGKr = p.PGK_Vm .* (rBPG .* rADP - rP3G .* rATP ./ p.PGK_Keq);
                v_PGKf = p.PGK_Vm .* (fBPG .* fADP - fP3G .* fATP ./ p.PGK_Keq);
                vObs = [-v_GAPDHr, v_GAPDHf];
%                 v_GAPDH = (-(p.TDH1_Vmr .* BPG .* NADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* NAD .* GAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./ p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
%                 v_PGK = p.PGK_Vm .* (BPG .* ADP - P3G .* ATP ./ p.PGK_Keq);
%                 vObs = -v_GAPDH;
            case 'gapdhr_s_revMM'
                % old-style
                p.TDH1_Keq=10.^x_temp(1).*0.0056; % []
                p.TDH1_Kgap=10.^x_temp(2).*2.48; % mM
                p.TDH1_Kbpg=10.^x_temp(3).*1.18; % mM
                p.TDH1_Knad=10.^x_temp(4).*2.92; %mM
                p.TDH1_Knadh=10.^x_temp(5).*0.022; % mM
                p.TDH1_Vm=10.^x_temp(6).*data.chosenVmax; % mM s^{-1}
                % check only one vmax
                if constantVm == 1
                    p.TDH1_Vm = data.chosenVmax; % fixing at the maximum
                else
                end
                    % % Standard method
                    % p.TDH1_Keq = x_temp(1); % []
                    % p.TDH1_Kgap = x_temp(2); % mM
                    % p.TDH1_Kbpg = x_temp(3); % mM
                    % p.TDH1_Knad = x_temp(4); %mM
                    % p.TDH1_Knadh = x_temp(5); % mM
                    % p.TDH1_Vm = x_temp(6).*data.chosenVmax; % mM s^{-1}

                    % calulate v (rateEquations)
                    v_GAPDHr = p.TDH1_Vm .* ( (BPG .* NADH) - (GAP .* NAD) ./ p.TDH1_Keq ) ./ (p.TDH1_Kbpg .* p.TDH1_Knadh) ./ ( (1 + GAP./p.TDH1_Kgap + BPG./p.TDH1_Kbpg).*(1 + NAD./p.TDH1_Knad +NADH./p.TDH1_Knadh) );
                    vObs = v_GAPDHr;
            otherwise
                disp('No kinetics have been selected');
        end
    otherwise
        disp('no enzyme selected in calcRates')
end
vextra = 1;
end