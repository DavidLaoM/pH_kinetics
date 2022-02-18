function [vObs, vextra] = calcRates(xtemp,simResult,data,setup)
enzymeName = setup.enzymeName;
constantVm = setup.constantVm;
ode = setup.ode;
sourceVm = setup.sourceVm;

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
                p.TDH1_Kgap = 10.^xtemp(2).*2.48; % mM
                p.TDH1_Kbpg = 10.^xtemp(3).*1.18; % mM
                p.TDH1_Knad = 10.^xtemp(4).*2.92; %mM
                p.TDH1_Knadh = 10.^xtemp(5).*0.022; % mM
                p.TDH1_Keq = data.chosenKeqGAPDH; % []
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
                % PGK
                p.PGK_Keq = data.chosenKeqPGK; % [] %/10
                p.PGK_Vm = 1306.45 / 60 * setup.excessPGK / data.chosenDF; % mM s^{-1} % corrected to make it appear in excess
                p.PGK_Katp = 0.3; % mM
                p.PGK_Kp3g = 0.53; % mM
                p.PGK_Kbpg = 0.003; % mM
                p.PGK_Kadp = 0.2; % mM
                % calulate v (rateEquations)
                v_GAPDH = (-(p.TDH1_Vmr .* BPG .* NADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* NAD .* GAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./ p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
                v_PGK = p.PGK_Vm .* (BPG .* ADP - P3G .* ATP ./ p.PGK_Keq);
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
            case 'gapdh_rev_simplified'
% % % %                 % overly simplified 2020-09-30
% % % %                 % recall parameters
% % % %                 p.TDH1_Vmr = data.chosenVmax * 10 .^ xtemp(1); %(UNIT!) 
% % % %                 p.GAPDH_Keq = data.chosenKeqGAPDH;
% % % %                 p.PGK_Vm = 1306.45 / 60;
% % % %                 p.PGK_Keq = data.chosenKeqPGK;                
% % % %                 % calculate the rates
% % % %                 v_GAPDHrev = p.TDH1_Vmr .* (BPG .* NADH ./ p.GAPDH_Keq - GAP .* NAD);
% % % %                 vObs = abs(v_GAPDHrev);
                % still keeping kms
                % recall parameters
                p.TDH1_Vmr = data.chosenVmax * 10 .^ xtemp(5); %(UNIT!) 
                p.GAPDH_Keq = data.chosenKeqGAPDH;
                p.PGK_Vm = 1306.45 / 60;
                p.PGK_Keq = data.chosenKeqPGK;      
                
                p.TDH1_Kgap=10.^xtemp(1).*2.48; % mM
                p.TDH1_Kbpg=10.^xtemp(2).*1.18; % mM
                p.TDH1_Knad=10.^xtemp(3).*2.92; %mM
                p.TDH1_Knadh=10.^xtemp(4).*0.022; % mM
                
%                 % calculate the rates
%                 v_GAPDHrev = ((p.TDH1_Vmr .* (BPG .* NADH ./ p.GAPDH_Keq - GAP .* NAD))./(p.TDH1_Kbpg .* p.TDH1_Knadh))./...
%                     ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./...
%                     p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
%                 vObs = abs(v_GAPDHrev);         
                
                % calculate the rates (adjusted)
                v_GAPDHrev = ((p.TDH1_Vmr .* (BPG .* NADH - GAP .* NAD  .* p.GAPDH_Keq))./(p.TDH1_Kbpg .* p.TDH1_Knadh))./...
                    ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./...
                    p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
                vObs = abs(v_GAPDHrev);       
% % % %                 % most simplified case
% % % %                 v_GAPDHrev = p.TDH1_Vmr * ones(size(v_GAPDHrev));
% % % %                 vObs = abs(v_GAPDHrev);
                
            otherwise
                disp('No kinetics have been selected');
        end
        
    case 'gapdh'
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
            case 'gapdh_fwd_simplified'
                % still keeping kms
                % recall parameters
                p.TDH1_Vmf = data.chosenVmax * 10 .^ xtemp(5); %(UNIT!) 
                p.GAPDH_Keq = data.chosenKeqGAPDH;
                p.PGK_Vm = 1306.45 / 60;
                p.PGK_Keq = data.chosenKeqPGK;      
                
                p.TDH1_Kgap=10.^xtemp(1).*2.48; % mM
                p.TDH1_Kbpg=10.^xtemp(2).*1.18; % mM
                p.TDH1_Knad=10.^xtemp(3).*2.92; %mM
                p.TDH1_Knadh=10.^xtemp(4).*0.022; % mM
                
                % calculate the rates
                v_GAPDHfwd = (p.TDH1_Vmf .* (GAP .* NAD - BPG .* NADH ./ p.GAPDH_Keq)./(p.TDH1_Kgap .* p.TDH1_Knad))./...
                    ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./...
                    p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
% % % %                 v_GAPDHfwd = (p.TDH1_Vmf .* (GAP .* NAD - BPG .* NADH ./ p.GAPDH_Keq));
                vObs = abs(v_GAPDHfwd);                
                
%                 % most simplified case
%                 v_GAPDHfwd = p.TDH1_Vmf * ones(size(v_GAPDHfwd));
%                 vObs = abs(v_GAPDHfwd);
                
            otherwise
                disp('No kinetics have been selected');
        end
        
    case 'eno'
        % metabolites
        P2G = simResult.y(:,1);
        PEP = simResult.y(:,2);
        % parameters
        p.ENO1_K2pg = 0.043 * 10.^xtemp(1); % mM
        p.ENO1_Kpep = 0.5 * 10.^xtemp(2); % mM
        p.ENO1_vm = data.chosenVmax * 10.^xtemp(3); % mM/s
        p.ENO1_Keq = data.chosenKeq; %[]
        % reactions
        v_ENO=((p.ENO1_vm./p.ENO1_K2pg).*(P2G-PEP./p.ENO1_Keq))./...
            (1+P2G./p.ENO1_K2pg+PEP./p.ENO1_Kpep);
        vObs = v_ENO; 
    case 'hxk'
        % metabolites
        GLCi = simResult.y(:,1); % mM
        G6P = simResult.y(:,2); % mM
        PG6 = simResult.y(:,3); % mM
        ATP = simResult.y(:,4); % mM
        ADP = simResult.y(:,5); % mM
        NADP = simResult.y(:,6); % mM
        NADPH = simResult.y(:,7); % mM
        % parameters
        p.HXK1_Kadp = 0.23 * 10 .^ xtemp(1); %mM
        p.HXK1_Katp = 0.15 * 10 .^ xtemp(2); %mM
        p.HXK1_Kg6p = 30 * 10 .^ xtemp(3); %mM
        p.HXK1_Kglc = 0.08 * 10 .^ xtemp(4); %mM
%         p.HXK1_Vm = .2*4.75* 10 .^ xtemp(5); %(UNIT!) 
        p.HXK1_Vm = data.chosenVmax * 10 .^ xtemp(5); %(UNIT!)
        p.HXK1_Kt6p = 0.2; %mM
        p.HXK1_Keq = data.chosenKeq_HXK; %[]
        p.G6PDH_Keq = data.chosenKeq_G6PDH; %[]
        p.G6PDH_Vm = 1306.45 / 60; % just copied from PGK since in that case it was big enough already
        % reactions
        v_GLK = (p.HXK1_Vm .* (ATP.*GLCi-((ADP.*G6P)./p.HXK1_Keq))) ./ ...
            ((p.HXK1_Katp.*p.HXK1_Kglc) .* (1+ATP./p.HXK1_Katp+ADP./p.HXK1_Kadp) .* ...
            (1+GLCi./p.HXK1_Kglc+G6P./p.HXK1_Kg6p+0./p.HXK1_Kt6p));
        v_G6PDH = p.G6PDH_Vm .* (NADP .* G6P - NADPH .* PG6 ./ p.G6PDH_Keq);
% % % %         vObs = v_G6PDH;
        vObs = cell(1,2);
        vObs{1} = v_GLK;
        vObs{2} = v_G6PDH;
    case 'ald'
        % metabolites
        FBP = simResult.y(:,1); % mM
        DHAP = simResult.y(:,2); % mM
        GAP = simResult.y(:,3); % mM
        G3P = simResult.y(:,4); % mM
        NADH = simResult.y(:,5); % mM
        NAD = simResult.y(:,6); % mM
        % parameters
        p.FBA1_Kf16bp = 0.451 * 10 .^ xtemp(1); % mM
        p.FBA1_Kglyceral3p = 2 * 10 .^ xtemp(2); % mM
        p.FBA1_Kdhap = 2.4 * 10 .^ xtemp(3); %mM
        p.FBA1_Vm = data.chosenVmax * 10 .^ xtemp(4); %(UNIT!) 
        p.FBA1_Keq = data.chosenKeq_FBA;
        p.TPI1_Vm = 1306.45 / 60;
        p.TPI1_Keq = data.chosenKeq_TPI;
        p.GPD1_Vm = 1306.45 / 60;
        p.GPD1_Keq = data.chosenKeq_GPD;
        % reactions
        v_ALD = p.FBA1_Vm .* (FBP-(GAP.*DHAP)./p.FBA1_Keq) ./...
            (p.FBA1_Kf16bp .* (FBP./p.FBA1_Kf16bp + (1 + GAP./p.FBA1_Kglyceral3p).*(1 + DHAP./p.FBA1_Kdhap)));
        v_GPD = p.GPD1_Vm .* (DHAP .* NADH - G3P .* NAD ./ p.GPD1_Keq);
        v_TPI = p.TPI1_Vm .* (DHAP - GAP ./ p.TPI1_Keq);
        vObs = [v_ALD, v_GPD, v_TPI];
    case 'pyk'
        % metabolites
        ADP = simResult.y(:,1); % mM
        NADH = simResult.y(:,2); % mM
        FBP = simResult.y(:,3); % mM
        PEP = simResult.y(:,4); % mM
        ATP = simResult.y(:,5); % mM
        NAD = simResult.y(:,6); % mM
        PYR = simResult.y(:,7); % mM
        LAC = simResult.y(:,8); % mM
        % parameters
        p.PYK1_Kadp = 10.^xtemp(1)  .*0.243; % mM
        p.PYK1_Katp = 10.^xtemp(2)  .*9.3; % mM
        p.PYK1_Kf16bp = 10.^xtemp(3).*0.2; % mM
        p.PYK1_Kpep = 10.^xtemp(4)  .*0.281; % mM
        p.PYK1_L = 10.^xtemp(5)     .*60000; % []
        p.PYK1_hill = 10.^xtemp(6)  .*4; % []
        p.PYK1_Vm = 10.^xtemp(7)    .*9.3167; % mM/s
        p.PYK1_Keq = data.chosenKeq_PYK;
        p.LDH1_Vm = 1306.45 / 60;
        p.LDH1_Keq = data.chosenKeq_LDH;
        % reactions
        v_PYK=(((p.PYK1_Vm./...
            (p.PYK1_Kadp.*p.PYK1_Kpep).*ADP.*PEP)./((1+ADP./p.PYK1_Kadp).*...
            (1+PEP./p.PYK1_Kpep))).*((PEP./p.PYK1_Kpep+1).^p.PYK1_hill./...
            (p.PYK1_L.*((ATP./p.PYK1_Katp+1)./(FBP./p.PYK1_Kf16bp+1)).^...
            p.PYK1_hill+(PEP./p.PYK1_Kpep+1).^p.PYK1_hill)));
        v_LDH = p.LDH1_Vm .*(PYR .* NADH - LAC .* NAD ./ p.LDH1_Keq);
        vObs = [v_PYK, v_LDH];
    case 'pgi'
        % metabolites
        G6P = simResult.y(:,1); % mM
        F6P = simResult.y(:,2); % mM
        ATP = simResult.y(:,3); % mM
        FBP = simResult.y(:,4); % mM
        ADP = simResult.y(:,5); % mM
        DHAP = simResult.y(:,6); % mM
        GAP = simResult.y(:,7); % mM
        NADH = simResult.y(:,8); % mM
        G3P = simResult.y(:,9); % mM
        NAD = simResult.y(:,10); % mM
        % parameters
        p.PGI1_Kg6p = 10 .^ xtemp(1) .* 1.0257; % mM
        p.PGI1_Kf6p = 10 .^ xtemp(2) .* 0.307; %mM
        p.PGI1_Vm = 10 .^ xtemp(3) .* data.chosenVmax;%[mM s^{-1}]
        p.PGI1_Keq = data.chosenKeq_PGI;
        p.PFK1_Vm = 1306.45 / 60;
        p.PFK1_Keq = data.chosenKeq_PFK;
        p.FBA1_Vm = 1306.45 / 60;
        p.FBA1_Keq = data.chosenKeq_FBA;
        p.TPI1_Vm = 1306.45 / 60;
        p.TPI1_Keq = data.chosenKeq_TPI;
        p.GPD1_Vm = 1306.45 / 60;
%         p.GPD1_Vm = 10 .^ xtemp(4) .* 1306.45 / 60;
        p.GPD1_Keq = data.chosenKeq_GPD;
        % reactions
        v_PGI = ((p.PGI1_Vm).*(G6P-(F6P./p.PGI1_Keq)))./...
            (1+G6P./p.PGI1_Kg6p+1+F6P./p.PGI1_Kf6p-1);
        v_PFK = p.PFK1_Vm .* (F6P .* ATP - FBP .* ADP ./ p.PGI1_Keq);
        v_ALD = p.FBA1_Vm .* (FBP - (GAP .* DHAP) ./ p.FBA1_Keq);
        v_TPI = p.TPI1_Vm .* (DHAP - GAP ./ p.TPI1_Keq);
        v_GPD = p.GPD1_Vm .* (DHAP .* NADH - G3P .* NAD ./ p.GPD1_Keq);
        vObs = [v_PGI, v_PFK, v_ALD, v_TPI, v_GPD];
    case 'pdc'
        % metabolites
        PYR = simResult.y(:,1); % mM
        CO2 = simResult.y(:,2); % mM
        AcAld = simResult.y(:,3); % mM
        ETOH = simResult.y(:,4); % mM
        NADH = simResult.y(:,5); % mM
        NAD = simResult.y(:,6); % mM
        % parameters
        p.PDC1_Kpyr = 10 .^ xtemp(1) .* 8.5; % mM
        p.PDC1_hill = 10 .^ xtemp(2) .* 1.9; % []
        p.PDC1_vm = 10.^xtemp(3) .* data.chosenVmax; % []
        p.ADH_Vm = 1306.45 / 60; % / 1000;
        p.ADH_Keq = data.chosenKeq_ADH;
        % reactions
        v_PDC = (p.PDC1_vm.*(PYR./p.PDC1_Kpyr).^p.PDC1_hill)./...
                (1+(PYR./p.PDC1_Kpyr).^p.PDC1_hill);
        v_ADH = p.ADH_Vm .* (AcAld .* NADH - ETOH .* NAD ./ p.ADH_Keq);
        vObs = [v_PDC, v_ADH];
    case 'pgm'
        % metabolites
        P3G = simResult.y(:,1); % mM
        P2G = simResult.y(:,2); % mM
        PEP = simResult.y(:,3); % mM
        PYR = simResult.y(:,4); % mM
        LAC = simResult.y(:,5); % mM
        ADP = simResult.y(:,6); % mM
        ATP = simResult.y(:,7); % mM
        NADH = simResult.y(:,8); % mM
        NAD = simResult.y(:,9); % mM
        % parameters
        p.GPM1_K2pg = 10 .^ xtemp(1) .* 0.08;
        p.GPM1_K3pg = 10 .^ xtemp(2) .* 1.2;
        p.GPM1_vm = 10.^xtemp(3) .* data.chosenVmax; % []
        p.GPM1_Keq = data.chosenKeq_PGM;
        p.ENO_Vm = 1306.45 / 60; % / 1000;
        p.ENO_Keq = data.chosenKeq_ENO;
        p.PYK_Vm = 1306.45 / 60; % / 1000;
        p.PYK_Keq = data.chosenKeq_PYK;
        p.LDH1_Vm = 1306.45 / 60; % / 1000;
        p.LDH1_Keq = data.chosenKeq_LDH;
        % reactions
        v_PGM = ((p.GPM1_vm ./ p.GPM1_K3pg) .* (P3G - P2G ./ p.GPM1_Keq))./...
            (1 + P3G ./ p.GPM1_K3pg + P2G ./ p.GPM1_K2pg);
        v_ENO = p.ENO_Vm .* (P2G - PEP ./ p.ENO_Keq);
        v_PYK = p.PYK_Vm .* (PEP .* ADP - PYR .* ATP ./ p.PYK_Keq);
        v_LDH = p.LDH1_Vm .*(PYR .* NADH - LAC .* NAD ./ p.LDH1_Keq);
        vObs = [v_PGM, v_ENO, v_PYK, v_LDH];
    case 'tpi'
        % metabolites
        FBP = simResult.y(:,1); % mM % just write down zero there
        DHAP = simResult.y(:,2); % mM
        GAP = simResult.y(:,3); % mM
        G3P = simResult.y(:,4); % mM
        NADH = simResult.y(:,5); % mM
        NAD = simResult.y(:,6); % mM
        % parameters
        p.TPI1_Kdhap = 10 .^ xtemp(1) .* 6.45; % mM
        p.TPI1_Kglyceral3p = 10 .^ xtemp(2) .* 5.25; % mM
        p.TPI1_Vm = 10.^xtemp(3) .* data.chosenVmax; % []
        p.TPI1_Keq = data.chosenKeq_TPI;
        p.GPD1_Vm = 1306.45 / 60; % / 1000;
        p.GPD1_Keq = data.chosenKeq_GPD;
        % reactions
        v_GPD = p.GPD1_Vm .* (DHAP .* NADH - G3P .* NAD ./ p.GPD1_Keq);
        v_TPI = (p.TPI1_Vm./p.TPI1_Kdhap.*(DHAP-GAP./p.TPI1_Keq))./...
            (1+DHAP./p.TPI1_Kdhap+GAP./p.TPI1_Kglyceral3p);
        v_ALD = zeros(size(v_GPD));
        vObs = [v_TPI, v_GPD, v_ALD];
    case 'pfk'
        % metabolites
        FBP = simResult.y(:,1); F16BP = FBP;% mM
        DHAP = simResult.y(:,2); % mM
        GAP = simResult.y(:,3); % mM
        G3P = simResult.y(:,4); % mM
        NADH = simResult.y(:,5); % mM
        NAD = simResult.y(:,6); % mM
        F6P = simResult.y(:,7); % mM
        ATP = simResult.y(:,8); % mM
        ADP = simResult.y(:,9); % mM
        F26BP = 0.1; % mM
        AMP = 0;
        % parameters
%         p.PFK_gR = 10 .^ xtemp(1) .* 5.12; %(1)
%         p.PFK_Kf6p = 10 .^ xtemp(2) .* 0.1; %(mmol/L)  %changed from 0.229 (Manchester kinetics)
%         p.PFK_Katp = 10 .^ xtemp(3) .* 0.71; %0.71; %(mmol/L) %changed from 0.068 (Manchester kinetics)
%         p.PFK_L = 10 .^ xtemp(4) .* 0.66; %(UNIT!)
%         p.PFK_Ciatp = 10 .^ xtemp(5) .* 100; %(1)
%         p.PFK_Kiatp = 10 .^ xtemp(6) .* 0.65; %(mmol/L)
%         p.PFK_Camp = 10 .^ xtemp(7) .* 0.0845; %(1) %CiAMP
%         p.PFK_Kamp = 10 .^ xtemp(8) .* 0.095; %(mmol/L)  %corrected value, new value from Teusink model
%         p.PFK_Cf26bp = 10 .^ xtemp(9) .* 0.0174; %(1) CiF26BP
%         p.PFK_Kf26bp = 10 .^ xtemp(10) .* 0.000682; %(mmol/L) %corrected value, new value from Teusink model
%         p.PFK_Cf16bp = 10 .^ xtemp(11) .* 0.397; %(1) CiF16BP
%         p.PFK_Kf16bp = 10 .^ xtemp(12) .* 0.111; %(mmol/L) %corrected value, new value from Teusink model
%         p.PFK_Catp = 10 .^ xtemp(13) .* 3; %(1)
            p.PFK.F26BP = 0.001; %adjusted from 0.02 (Teusink model), new value based on Canelas data
%         p.PFK_Vm = 10.^xtemp(14) .* data.chosenVmax; % []
        p.FBA1_Keq = data.chosenKeq_FBA;
        p.FBA1_Vm = 1306.45 / 60; % / 1000;
        p.GPD1_Keq = data.chosenKeq_GPD;
        p.GPD1_Vm = 1306.45 / 60; % / 1000;
        p.TPI1_Keq = data.chosenKeq_TPI;
        p.TPI1_Vm = 1306.45 / 60; % / 1000;
        p.PFK_Keq = data.chosenKeq_PFK;
        p.PFK_Vm = 10.^xtemp(1) .* data.chosenVmax; % []
        % reactions
%         for casePFK = 1
%             PFK_nom=(p.PFK_Vm.*p.PFK_gR.*(F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp).*(1+(F6P./p.PFK_Kf6p)+(ATP./p.PFK_Katp)+p.PFK_gR.*((F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp))));
%             PFK_denom=(1+F6P./p.PFK_Kf6p+ATP./p.PFK_Katp+(p.PFK_gR.*(F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp))).^2+...
%                 p.PFK_L.*...
%                 ((1+p.PFK_Ciatp.*(ATP./p.PFK_Kiatp))./(1+ATP./p.PFK_Kiatp)).^2.*...
%                 ((1+p.PFK_Camp.*(AMP./p.PFK_Kamp))./(1+AMP./p.PFK_Kamp)).^2.*...
%                 ((1+((p.PFK_Cf26bp*F26BP)./(p.PFK_Kf26bp))+((p.PFK_Cf16bp.*F16BP)./(p.PFK_Kf16bp)))./(1+(F26BP./p.PFK_Kf26bp)+(F16BP./p.PFK_Kf16bp))).^2.*...
%                 (1+p.PFK_Catp.*(ATP./p.PFK_Katp)).^2;
%             v_PFK = PFK_nom ./ PFK_denom;
%         end
        for casePFK = 2
            v_PFK = p.PFK_Vm .* (F6P .* ATP - (FBP .* ADP)./p.PFK_Keq);
        end
        v_ALD = p.FBA1_Vm .* (FBP-(GAP.*DHAP)./p.FBA1_Keq);
        v_GPD = p.GPD1_Vm .* (DHAP .* NADH - G3P .* NAD ./ p.GPD1_Keq);
        v_TPI = p.TPI1_Vm .* (DHAP - GAP ./ p.TPI1_Keq);
        vObs = [v_PFK, v_ALD, v_GPD, v_TPI];
    otherwise
        disp('no enzyme selected in calcRates')
end
vextra = 1;
end