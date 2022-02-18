function [simResult] = simSys(xtemp,data,setup)
% file where simulation is run
enzymeName = setup.enzymeName;
constantVm = setup.constantVm;
ode = setup.ode;
plotEachSim = setup.plotEachSim;
legendamets = setup.PSAmets;
sourceVm = setup.sourceVm;
% disp(data.i);        
switch enzymeName
    case 'gapdhr'
        switch ode
            case 'vanHeerden2014'
                odefun = @odeGAPDH_vHeerden;
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
        %         p.PGK_Vm = 1306.45 / 60 * setup.excessPGK; % / data.chosenDF; % mM s^{-1} % corrected to make it appear in excess
                p.PGK_Vm = 1306.45 / 60 * setup.excessPGK / data.chosenDF; % mM s^{-1} % corrected to make it appear in excess
                p.PGK_Katp = 0.3; % mM
                p.PGK_Kp3g = 0.53; % mM
                p.PGK_Kbpg = 0.003; % mM
                p.PGK_Kadp = 0.2; % mM
        %         %display check
        %         disp(p.TDH1_Vmf);
        %         disp(p.TDH1_Vmr);
        %         disp(p.PGK_Vm);
            case 'gapdhr_s_revMM'
                odefun = @odeGAPDHr;
                p.TDH1_Keq=10.^xtemp(1).*0.0056; % []
        % % % %         p.TDH1_Keq=setup.chosenKeqgapdh;
                p.TDH1_Kgap=10.^xtemp(2).*2.48; % mM
                p.TDH1_Kbpg=10.^xtemp(3).*1.18; % mM
                p.TDH1_Knad=10.^xtemp(4).*2.92; %mM
                p.TDH1_Knadh=10.^xtemp(5).*0.022; % mM
                p.TDH1_Vm=10.^xtemp(6).*data.chosenVmax; % mM s^{-1}
                p.linkreaction = 10.^xtemp(7)./data.chosenLink;
            case 'gapdh_rev_simplified'
% % % %                 % overly simplified 2020-09-30
% % % %                 odefun = @odeGAPDHrev_simplified;
% % % %                 p.TDH1_Vmr = data.chosenVmax * 10 .^ xtemp(1); %(UNIT!) 
% % % %                 p.GAPDH_Keq = data.chosenKeqGAPDH;
% % % %                 p.PGK_Vm = 1306.45 / 60;
% % % %                 p.PGK_Keq = data.chosenKeqPGK;
                % still keeping kms
                odefun = @odeGAPDHrev_simplified;
                p.TDH1_Kgap=10.^xtemp(1).*2.48; % mM
                p.TDH1_Kbpg=10.^xtemp(2).*1.18; % mM
                p.TDH1_Knad=10.^xtemp(3).*2.92; %mM
                p.TDH1_Knadh=10.^xtemp(4).*0.022; % mM
                p.TDH1_Vmr = data.chosenVmax * 10 .^ xtemp(5); %(UNIT!) 
                p.GAPDH_Keq = data.chosenKeqGAPDH;
                p.PGK_Vm = 1306.45 / 60;
%                 p.PGK_Vm = 1306.45 * 1306.45 / 60;
                p.PGK_Keq = data.chosenKeqPGK;
            otherwise
                disp('No ode being selected');
        end

        % check only one vmax
        if constantVm == 1
            p.TDH1_Vm = data.chosenVmax; % fixing at the maximum
        else
        end
        
        % set initial conditions and timespan
        P3Go = 5;
        ATPo = 1;
        BPGo = 0;
        ADPo = 0;
        NADo = 0;
        GAPo = 0;
        PHOSo = 500;
        NADHo = data.chosenNADini;
        
%         %testing in the PSA
%         P3Go = 4.5;
        
        tspan   = [0 300];  % time [300 s (6 min)]
% % % %         tspan = [0 3000];
        y0 = [P3Go ATPo BPGo ADPo NADo GAPo PHOSo NADHo];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
        % run ode15s
        f = 1; % not being used for PSA
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);

        simResult.t = t;
        simResult.y = y;
        
    case 'gapdh'
        switch ode
            case 'gapdh_fwd_simplified'
                % still keeping kms
                odefun = @odeGAPDHfwd_simplified;
                p.TDH1_Kgap=10.^xtemp(1).*2.48; % mM
                p.TDH1_Kbpg=10.^xtemp(2).*1.18; % mM
                p.TDH1_Knad=10.^xtemp(3).*2.92; %mM
                p.TDH1_Knadh=10.^xtemp(4).*0.022; % mM
                p.TDH1_Vmf = data.chosenVmax * 10 .^ xtemp(5); %(UNIT!) 
                p.GAPDH_Keq = data.chosenKeqGAPDH;
                p.PGK_Vm = 1306.45 / 60;
                p.PGK_Keq = data.chosenKeqPGK;
            otherwise
                disp('No ode being selected');
        end

        % check only one vmax
        if constantVm == 1
            p.TDH1_Vm = data.chosenVmax; % fixing at the maximum
        else
        end
        
        % set initial conditions and timespan
        P3Go = 0;
        ATPo = 0;
        BPGo = 0;
        ADPo = 10;
        NADo = 1;
        GAPo = 5.8; %0.12;
        PHOSo = 0;
        NADHo = data.chosenNADini;

        tspan   = [0 600];  % time [300 s (6 min)]
% % % %         tspan = [0 3000];
        y0 = [P3Go ATPo BPGo ADPo NADo GAPo PHOSo NADHo];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
        % run ode15s
        f = 1; % not being used for PSA
        % % 
        % addition for DoE
        % DLM 2020-11-19
        if isfield(setup,'PSAstudy') % PSA?
            % which PSA?
            if isfield(setup,'PSAstudy_GAPDH_kgap')
                if setup.PSAstudy_GAPDH_kgap == 1
                    y0(6) = setup.PSAvals.GAPDH_kgap_ini(setup.idx);
                    tspan = [0 1000];
%                     y0(1) = 0.6;
                end
            end
        end
        % % 
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);

        simResult.t = t;
        simResult.y = y;
        
    case 'eno'
        % setParameterStructure
        p.ENO1_K2pg = 0.043 * 10.^xtemp(1); % mM
        p.ENO1_Kpep = 0.5 * 10.^xtemp(2); % mM
        p.ENO1_vm = data.chosenVmax * 10.^xtemp(3); % mM/s
        p.ENO1_Keq = data.chosenKeq; %[]
        % simulate
        P2Go = 6; % reference from 'in_vitro_reactions.pdf'
        PEPo = data.chosenPEPini;
        tspan = [0 300];
        y0 = [P2Go PEPo];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));%odefun = @odeGAPDHr;
        f = 1; % not being used for PSA
        odefun = @odeENO;
        % % 
        % addition for DoE
        % DLM 2020-11-19
        if isfield(setup,'PSAstudy') % PSA?
            % which PSA?
            if isfield(setup,'PSAstudy_ENO_k2pg')
                if setup.PSAstudy_ENO_k2pg == 1
                    y0(1) = setup.PSAvals.ENO_k2pg_ini(setup.idx);
                    tspan = [0 1000];
%                     y0(1) = 0.6;
                end
            end
            if isfield(setup,'PSAstudy_ENO_kpep')
                if setup.PSAstudy_ENO_kpep == 1
                    tspan = [0 1000];
% % % %                     y0(1) = 0.6;
                end
            end
            y0(2) = 0;
        end
        % % 
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);  
        simResult.t = t;
        simResult.y = y;
    case 'hxk'
        % setParameterStructure
        p.HXK1_Kadp = 0.23 * 10 .^ xtemp(1); %mM
        p.HXK1_Katp = 0.15 * 10 .^ xtemp(2); %mM
        p.HXK1_Kg6p = 30 * 10 .^ xtemp(3); %mM
        p.HXK1_Kglc = 0.08 * 10 .^ xtemp(4); %mM
        p.HXK1_Vm = data.chosenVmax * 10 .^ xtemp(5); %(UNIT!) 
        p.HXK1_Kt6p = 0.2; %mM
        p.HXK1_Keq = data.chosenKeq_HXK; %[]
        p.G6PDH_Keq = data.chosenKeq_G6PDH; %[]
        p.G6PDH_Vm = 1306.45 / 60; %/ 100000; % just copied from PGK since in that case it was big enough already
% % % %         p.G6PDH_Vm = 130600000.45 / 60; %/ 100000; % just copied from PGK since in that case it was big enough already
        % simulate
        Glco = 10; % mM
        G6Po = 0; % mM
        PG6o = 0; % mM
        ATPo = 1; % mM
        ADPo = 0; % mM
        NADPo = 1; % mM
        NADPHo = data.chosenNADPHini; % mM
        tspan = [0 300];
        y0 = [Glco G6Po PG6o ATPo ADPo NADPo NADPHo];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
        f = 1; % not being used for PSA
        odefun = @odeHXK;
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);  
        simResult.t = t;
        simResult.y = y;
    case 'ald'
        % setParameterStructure
        p.FBA1_Kf16bp = 0.451 * 10 .^ xtemp(1); % mM
        p.FBA1_Kglyceral3p = 2 * 10 .^ xtemp(2); % mM
        p.FBA1_Kdhap = 2.4 * 10 .^ xtemp(3); %mM
        p.FBA1_Vm = data.chosenVmax * 10 .^ xtemp(4); %(UNIT!) 
        p.FBA1_Keq = data.chosenKeq_FBA;
% % % %         p.TPI1_Vm = 1306.45 / 60 * 10 .^ xtemp(5);% /data.chosenDF;% / 10;
        p.TPI1_Vm = 1306.45 / 60;% /data.chosenDF;% / 10;
        p.TPI1_Keq = data.chosenKeq_TPI;
%         p.GPD1_Vm = 1306.45 / 60 * 10 .^ xtemp(6);% /data.chosenDF;% / 100;
% % % %         p.GPD1_Vm = 1306.45 / 60 / 50 * 10 .^ xtemp(6);% /data.chosenDF;% / 100;
        p.GPD1_Vm = 1306.45 / 60;% /data.chosenDF;% / 100;
        p.GPD1_Keq = data.chosenKeq_GPD;
        % simulate
        FBPo = 2; % mM
        DHAPo = 0; %mM
        GAPo = 0; %mM
        G3Po = 0; %mM
        NADHo = data.chosenNADHini; %mM
        NADo = 0; %mM
        tspan = [0 420];
        y0 = [FBPo DHAPo GAPo G3Po NADHo NADo];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
        f = 1; % not being used for PSA
        odefun = @odeALD;
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);  
        simResult.t = t;
        simResult.y = y;
    case 'pyk'
        % setParameterStructure
        p.PYK1_Kadp = 10.^xtemp(1)  .*0.243; % mM
        p.PYK1_Katp = 10.^xtemp(2)  .*9.3; % mM
        p.PYK1_Kf16bp = 10.^xtemp(3).*0.2; % mM
        p.PYK1_Kpep = 10.^xtemp(4)  .*0.281; % mM
        p.PYK1_L = 10.^xtemp(5)     .*60000; % []
        p.PYK1_hill = 10.^xtemp(6)  .*4; % []
        p.PYK1_Vm = 10.^xtemp(7)    .*data.chosenVmax; % mM/s
        p.PYK1_Keq = data.chosenKeq_PYK;
        p.LDH1_Vm = 1306.45 / 60; % / 1000;
        p.LDH1_Keq = data.chosenKeq_LDH;
        % simulate
        ADPo = 10;
        NADHo = data.chosenNADHini;
        FBPo = 1;
        PEPo = 2;
        ATPo = 0;
        NADo = 0;
        PYRo = 0;
        LACo = 0;
        tspan = [0 600];
        y0 = [ADPo NADHo FBPo PEPo ATPo NADo PYRo LACo];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
        f = 1; % not being used for PSA
% % % %         % % % %
% % % %         disp(p.PYK1_Vm); % % % %
% % % %         % % % %
        odefun = @odePYK;
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);  
        simResult.t = t;
        simResult.y = y;
    case 'pgi'
        % setParameterStructure
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
        % simulate
        G6Po = 5;
        F6Po = 0;
        ATPo = 1; 
        FBPo = 0;
        ADPo = 0;
        DHAPo = 0;
        GAPo = 0;
        NADHo = data.chosenNADHini;
        G3Po = 0;
        NADo = 0;
        tspan = [0 1000];
        y0 = [G6Po F6Po ATPo FBPo ADPo DHAPo GAPo NADHo G3Po NADo];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
        f = 1; % not being used for PSA
        odefun = @odePGI;
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);  
        simResult.t = t;
        simResult.y = y;
    case 'pdc'
        % setParameterStructure
        p.PDC1_Kpyr = 10 .^ xtemp(1) .* 8.5; % mM
        p.PDC1_hill = 10 .^ xtemp(2) .* 1.9; % []
        p.PDC1_vm = 10.^xtemp(3) .* data.chosenVmax; % []
        p.ADH_Vm = 1306.45 / 60; % / 1000;
        p.ADH_Keq = data.chosenKeq_ADH;
        % simulate
        PYRo = 50; % mM
        CO2o = 0; % mM
        AcAldo = 0; % mM
        ETOHo = 0; % mM
        NADHo = data.chosenNADHini; % mM
        NADo = 0; % mM
        tspan = [0 600];
        y0 = [PYRo CO2o AcAldo ETOHo NADHo NADo];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
        f = 1; % not being used for PSA
        odefun = @odePDC;
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);  
        simResult.t = t;
        simResult.y = y;
    case 'pgm'
        % setParameterStructure
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
        % simulate
        P3Go = 5; % mM
        P2Go = 0; % mM
        PEPo = 0; % mM
        PYRo = 0; % mM
        LACo = 0; % mM
        ADPo = 10; % mM
        ATPo = 0; % mM
        NADHo = data.chosenNADHini; % mM
        NADo = 0; % mM
        tspan = [0 600];
        y0 = [P3Go P2Go PEPo PYRo LACo ADPo ATPo NADHo NADo];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
        f = 1; % not being used for PSA
        odefun = @odePGM;
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);  
        simResult.t = t;
        simResult.y = y;
    case 'tpi'
        % setParameterStructure
        p.TPI1_Kdhap = 10 .^ xtemp(1) .* 6.45; % mM
        p.TPI1_Kglyceral3p = 10 .^ xtemp(2) .* 5.25; % mM
        p.TPI1_Vm = 10.^xtemp(3) .* data.chosenVmax; % []
        p.TPI1_Keq = data.chosenKeq_TPI;
        p.GPD1_Vm = 1306.45 / 60; % / 1000;
        p.GPD1_Keq = data.chosenKeq_GPD;
        % simulate
        FBPo = 0; % mM
        DHAPo = 0; % mM
        GAPo = 5.8; % mM
        G3Po = 0; % mM
        NADHo = data.chosenNADHini; % mM
        NADo = 0; % mM
        tspan = [0 600];
        y0 = [FBPo DHAPo GAPo G3Po NADHo NADo];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
        f = 1; % not being used for PSA
        odefun = @odeTPI;
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);  
        simResult.t = t;
        simResult.y = y;
    case 'pfk'
        % setParameterStructure
        
        % Problem reaction PFK: selecting the right parameters
        problemStudy = setup.problemStudy;
        switch problemStudy
            case 'onlyVmax'
                p.PFK.F26BP = 0.001; %adjusted from 0.02 (Teusink model), new value based on Canelas data
                p.PFK_Vm = 10.^xtemp(1) .* data.chosenVmax; % []
            case 'fullKinetics_paramsPartFixed'
                p.PFK_gR = 10 .^ xtemp(1) .* 5.12; %(1)
                p.PFK_Kf6p = 10 .^ xtemp(2) .* 0.1; %(mmol/L)  %changed from 0.229 (Manchester kinetics)
                p.PFK_Katp = 10 .^ xtemp(3) .* 0.71; %0.71; %(mmol/L) %changed from 0.068 (Manchester kinetics)
                p.PFK_L = 10 .^ xtemp(4) .* 0.66; %(UNIT!)
                p.PFK_Ciatp = 10 .^ xtemp(5) .* 100; %(1)

                p.PFK_Kiatp = 10 .^ xtemp(6) .* 0.65; %(mmol/L)
                p.PFK_Camp = 10 .^ xtemp(7) .* 0.0845; %(1) %CiAMP
                p.PFK_Kamp = 10 .^ xtemp(8) .* 0.095; %(mmol/L)  %corrected value, new value from Teusink model
                p.PFK_Cf26bp = 10 .^ xtemp(9) .* 0.0174; %(1) CiF26BP
                p.PFK_Kf26bp = 10 .^ xtemp(10) .* 0.000682; %(mmol/L) %corrected value, new value from Teusink model

                p.PFK_Cf16bp = 10 .^ xtemp(11) .* 0.397; %(1) CiF16BP
                p.PFK_Kf16bp = 10 .^ xtemp(12) .* 0.111; %(mmol/L) %corrected value, new value from Teusink model
                p.PFK_Catp = 10 .^ xtemp(13) .* 3; %(1)

                p.PFK.F26BP = 0.001; %adjusted from 0.02 (Teusink model), new value based on Canelas data

                p.PFK_Vm = 10.^xtemp(14) .* data.chosenVmax; % []
            case 'fullKinetics_paramsAllFlexible'
                p.PFK_gR = 10 .^ xtemp(1) .* 5.12; %(1)
                p.PFK_Kf6p = 10 .^ xtemp(2) .* 0.1; %(mmol/L)  %changed from 0.229 (Manchester kinetics)
                p.PFK_Katp = 10 .^ xtemp(3) .* 0.71; %0.71; %(mmol/L) %changed from 0.068 (Manchester kinetics)
                p.PFK_L = 10 .^ xtemp(4) .* 0.66; %(UNIT!)
                p.PFK_Ciatp = 10 .^ xtemp(5) .* 100; %(1)

                p.PFK_Kiatp = 10 .^ xtemp(6) .* 0.65; %(mmol/L)
                p.PFK_Camp = 10 .^ xtemp(7) .* 0.0845; %(1) %CiAMP
                p.PFK_Kamp = 10 .^ xtemp(8) .* 0.095; %(mmol/L)  %corrected value, new value from Teusink model
                p.PFK_Cf26bp = 10 .^ xtemp(9) .* 0.0174; %(1) CiF26BP
                p.PFK_Kf26bp = 10 .^ xtemp(10) .* 0.000682; %(mmol/L) %corrected value, new value from Teusink model

                p.PFK_Cf16bp = 10 .^ xtemp(11) .* 0.397; %(1) CiF16BP
                p.PFK_Kf16bp = 10 .^ xtemp(12) .* 0.111; %(mmol/L) %corrected value, new value from Teusink model
                p.PFK_Catp = 10 .^ xtemp(13) .* 3; %(1)

                p.PFK.F26BP = 0.001; %adjusted from 0.02 (Teusink model), new value based on Canelas data
                p.PFK_Vm = 10.^xtemp(14) .* data.chosenVmax; % []
                
             otherwise
                disp('Warning: problem study (reaction kinetics) have not been selected');
        end
        p.PFK_Keq = data.chosenKeq_PFK;
        
        % linking reactions
        p.FBA1_Keq = data.chosenKeq_FBA;
        p.FBA1_Vm = 1306.45 / 60; % / 1000;
        p.GPD1_Keq = data.chosenKeq_GPD;
        p.GPD1_Vm = 1306.45 / 60; % / 1000;
        p.TPI1_Keq = data.chosenKeq_TPI;
        p.TPI1_Vm = 1306.45 / 60; % / 1000;
        % simulate
        FBPo = 0; % mM
        DHAPo = 0; % mM
        GAPo = 0; % mM
        G3Po = 0; % mM
        NADHo = data.chosenNADHini; % mM
        NADo = 0; % mM
        F6Po = 10; % or 0.25? % mM
% % % %         F6Po = 0.25; % mM
        ATPo = 0.5; % mM
        ADPo = 0; % mM
%         tspan = [0 600];
%         tspan = [0 1000];
        tspan = [0 2000];
        y0 = [FBPo DHAPo GAPo G3Po NADHo NADo F6Po ATPo ADPo];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
        f = 1; % not being used for PSA
        odefun = @odePFK;        
        % % 
        % addition for DoE
        % DLM 2020-11-19
        if isfield(setup,'PSAstudy') % PSA?
            % which PSA?
            if isfield(setup,'PSAstudy_PFK_kf6p')
                if setup.PSAstudy_PFK_kf6p == 1
                    y0(7) = setup.PSAvals.PFK_kf6p_ini(setup.idx);
%                     tspan = [0 1000];
%                     y0(1) = 0.6;
                end
            end
        end
        % % 
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);  
        simResult.t = t;
        simResult.y = y;
    otherwise
        disp('no enzyme has been selected in simSys');
end

if plotEachSim == 2
    [~,n] = size(y);
    figure
    for i = 1:n
        subplot(3,3,i)
        plot(t,y(:,i),'b')
        title(legendamets{i})
        hold on
    end
    suptitle('Simple plot system simulations')
end 

end

% % Scheme PGI
% metabolites involved
% name	number	starting concentration(mM)	mass balance
% G6P	x1      5                           - v_PGI
% F6P	x2      0                           + v_PGI - v_PFK
% ATP	x3      1                           - v_PFK
% FBP	x4      0                           + v_PFK - v_ALD
% ADP	x5      0                           + v_PFK
% DHAP	x6      0                           + v_ALD - v_TPI - v_GPD
% GAP	x7      0                           + v_ALD + v_TPI
% NADH	x8      0.15                        - v_GPD
% G3P	x9      0                           + v_GPD
% NAD	x10     0                           + v_GPD
% 
% reaction kinetics
% v_PGI = ((p.PGI1_Vm).*(G6P-(F6P./p.PGI1_Keq)))./...
%     (1+G6P./p.PGI1_Kg6p+1+F6P./p.PGI1_Kf6p-1);
% v_PFK = p.PFK1_Vm .* (F6P .* ATP - FBP .* ADP ./ p.PGI1_Keq);
% v_ALD = p.FBA1_Vm .* (FBP - (GAP .* DHAP) ./ p.FBA1_Keq);
% v_TPI = p.TPI1_Vm .* (DHAP - GAP ./ p.TPI1_Keq);
% v_GPD = p.GPD1_Vm .* (DHAP .* NADH - G3P .* NAD ./ p.GPD1_Keq);
% 
% parameters of interest are
% name            number
% p.PGI1_Keq      1
% p.PGI1_Kg6p     2
% p.PGI1_Kf6p     3
% p.PGI1_Vm       4-15

% % memoryDump ('old' new kinetics. Wrong location for keqs.)
% 
% switch ode
%     case 'vanHeerden2014'
%         odefun = @odeGAPDH_vHeerden;
%         % GAPDH
% % % % %         p.TDH1_Vmf = 10.^xtemp(1).*data.chosenVmf;% mM s^{-1}
% % % % %         p.TDH1_Vmf = 10.^xtemp(1).*1184.52/60;% mM s^{-1}
%         p.TDH1_Vmf = 10.^xtemp(1).*1184.52/60 / data.chosenDF;% mM s^{-1}        
%         p.TDH1_Kgap = 10.^xtemp(2).*2.48; % mM
%         p.TDH1_Kbpg = 10.^xtemp(3).*1.18; % mM
%         p.TDH1_Knad = 10.^xtemp(4).*2.92; %mM
%         p.TDH1_Knadh = 10.^xtemp(5).*0.022; % mM
% % % % %         p.TDH1_Vmr = 10.^xtemp(6).*data.chosenVmr; % mM s^{-1} %.*data.chosenVmax
% % % % %         p.TDH1_Vmr = 10.^xtemp(6).*6549.8/60; % mM s^{-1} %.*data.chosenVmax
%         p.TDH1_Vmr = 10.^xtemp(6).*6549.8/60 / data.chosenDF; % mM s^{-1} %.*data.chosenVmax
%         p.TDH1_Keq = data.chosenKeqGAPDH; %=10.^xtemp(1).*0.0056; % []
% % % % %         p.TDH1_Keq = 0.0056; %=10.^xtemp(1).*0.0056; % []
% %         p.TDH1_Vmr = 10.^xtemp(6).*0.0021; % mM s^{-1} %.*data.chosenVmax       % % % %
% %         p.TDH1_Keq = 10.^xtemp(6).*0.0056; %=10.^xtemp(1).*0.0056; % []         % % % %
%         % PGK
%         p.PGK_Keq = data.chosenKeqPGK; % [] %/10
% % % % %         p.PGK_Keq = 3200; % [] %/10
% % % % %         p.PGK_Keq = 10.^xtemp(7) .* 3200; % []
%         p.PGK_Vm = 1306.45/60*1; % mM s^{-1}
%         p.PGK_Katp = 0.3; % mM
%         p.PGK_Kp3g = 0.53; % mM
%         p.PGK_Kbpg = 0.003; % mM
%         p.PGK_Kadp = 0.2; % mM
%     case 'gapdhr_s_revMM'
%         odefun = @odeGAPDHr;
%         p.TDH1_Keq=10.^xtemp(1).*0.0056; % []
%         p.TDH1_Kgap=10.^xtemp(2).*2.48; % mM
%         p.TDH1_Kbpg=10.^xtemp(3).*1.18; % mM
%         p.TDH1_Knad=10.^xtemp(4).*2.92; %mM
%         p.TDH1_Knadh=10.^xtemp(5).*0.022; % mM
%         p.TDH1_Vm=10.^xtemp(6).*data.chosenVmax; % mM s^{-1}
%         p.linkreaction = 10.^xtemp(7)./data.chosenLink;
%     otherwise
%         disp('No ode being selected');
% end


