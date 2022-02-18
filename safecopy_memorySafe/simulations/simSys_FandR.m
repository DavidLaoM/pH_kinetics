function [simResult] = simSys_FandR(xtemp,data,setup)
% file where simulation is run
enzymeName = setup.enzymeName;
constantVm = setup.constantVm;
ode = setup.ode;
plotEachSim = setup.plotEachSim;
legendamets = setup.PSAmets;
sourceVm = setup.sourceVm;
typeVm = setup.typeVm;
ode_pH = setup.ode_pH;

% disp(data.i);
switch ode
    case 'vanHeerden2014'
        odefun = @odeGAPDH_vHeerden_FandR;
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
                        switch ode_pH
                            case 'on_revMM_bitri'
                                p.TDH1_Vmr_R = 10.^xtemp(6).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                                p.TDH1_Vmf_R = 10.^xtemp(1).*30.98;
%                                 p.TDH1_Vmr_R = 10.^xtemp(6).*p.TDH1_Keq; % mM s^{-1}
%                                 
%                                 p.TDH1_Vmf_R
%                                 
%                                 p.TDH1_Vmr_R = 
                            otherwise
                                p.TDH1_Vmr_R = 10.^xtemp(6).*setup.exp_vmax_gapdhr(6) / data.chosenDF; % mM s^{-1}
                        end
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
    otherwise
        disp('No ode being selected');
end
        
switch enzymeName
    case 'gapdhr'
        % check only one vmax
        if constantVm == 1
            p.TDH1_Vm = data.chosenVmax; % fixing at the maximum
        else
        end
        
        % set initial conditions and timespan
        rP3Go = 5;
        rATPo = 1;
        rBPGo = 0;
        rADPo = 0;
        rNADo = 0;
        rGAPo = 0;
        rPHOSo = 500;
        rNADHo = data.gapdhR.chosenNADini;

        fP3Go = 0;
        fATPo = 0;
        fBPGo = 0;
        fADPo = 10;
        fNADo = 1;
        fGAPo = 5.8;
        fPHOSo = 500;
        fNADHo = data.gapdhF.chosenNADini;        
 
        y0r = [rP3Go rATPo rBPGo rADPo rNADo rGAPo rPHOSo rNADHo];
        y0f = [fP3Go fATPo fBPGo fADPo fNADo fGAPo fPHOSo fNADHo];
        
        y0 = [y0r, y0f];
        
        tspan = [0 600];
        options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
        % run ode15s
        f = 1; % not being used for PSA
        [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);

        simResult.t = t;
        simResult.y = y;
        
    otherwise
        disp('no enzyme has been selected in simSys');
end

if plotEachSim == 1
    [~,n] = size(y);
    n12 = n/2;
    figure
    for i = 1:n12
        subplot(3,3,i)
        plot(t,y(:,i),'b')
        title(legendamets{i})
        hold on
    end
    suptitle('System simulation: Reverse direction')
    figure
    for i = (n12+1):(n12*2)
        subplot(3,3,i-n12)
        plot(t,y(:,i),'b')
        title(legendamets{i-n12})
        hold on
    end
    suptitle('System simulation: Forward direction')
end 

end

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