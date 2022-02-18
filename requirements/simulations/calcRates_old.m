function [vObs, vextra] = calcRates(x_temp,simResult,data,setup)
enzymeName = setup.enzymeName;
constantVm = setup.constantVm;
ode = setup.ode;

switch enzymeName
    case 'gapdhr'
        % metabolites
        P3G = simResult.y(:,1);
        ATP = simResult.y(:,2);
        BPG = simResult.y(:,3);
        ADP = simResult.y(:,4);
        NAD = simResult.y(:,5);
        GAP = simResult.y(:,6);
        PHOS = simResult.y(:,7);
        NADH = simResult.y(:,8);
        
        switch ode
            case 'vanHeerden2014'
                % GAPDH
                p.TDH1_Vmf = 10.^x_temp(1).*data.chosenVmf;% mM s^{-1}
                p.TDH1_Kgap = 10.^x_temp(2).*2.48; % mM
                p.TDH1_Kbpg = 10.^x_temp(3).*1.18; % mM
                p.TDH1_Knad = 10.^x_temp(4).*2.92; %mM
                p.TDH1_Knadh = 10.^x_temp(5).*0.022; % mM
                p.TDH1_Vmr = 10.^x_temp(6).*data.chosenVmr; % mM s^{-1} %.*data.chosenVmax
                p.TDH1_Keq = data.chosenKeqGAPDH; %=10.^xtemp(1).*0.0056; % []
                % PGK
                p.PGK_Vm = 1306.45/60*1; % mM s^{-1}
                p.PGK_Keq = data.chosenKeqPGK; % [] %/10
                p.PGK_Katp = 0.3; % mM
                p.PGK_Kp3g = 0.53; % mM
                p.PGK_Kbpg = 0.003; % mM
                p.PGK_Kadp = 0.2; % mM
                % calulate v (rateEquations)
                v_GAPDH = (-(p.TDH1_Vmr .* BPG .* NADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* NAD .* PHOS .* GAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + PHOS) .* (1 + BPG ./ p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
                v_PGK = p.PGK_Vm .* (p.PGK_Keq .* BPG .* ADP - P3G .* ATP) ./ (p.PGK_Katp .* p.PGK_Kp3g .* (1 + BPG ./ p.PGK_Kbpg + P3G ./ p.PGK_Kp3g) .* (1 + ADP ./ p.PGK_Kadp + ATP ./ p.PGK_Katp));

                vObs = -v_GAPDH;
                
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