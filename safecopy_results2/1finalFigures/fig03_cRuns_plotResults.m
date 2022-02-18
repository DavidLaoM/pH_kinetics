%1 Get data
%2 Make array together
%3 Plot
%4 Save array
%5 Save plot

%% 1 get data
% vf
load('output_vf_run1.mat'); output_vf_run1 = output; clear output
load('output_vf_run2.mat'); output_vf_run2 = output; clear output
load('output_vf_run3.mat'); output_vf_run3 = output; clear output
load('output_vf_run4.mat'); output_vf_run4 = output; clear output
load('output_vf_run5.mat'); output_vf_run5 = output; clear output
load('output_vf_run6.mat'); output_vf_run6 = output; clear output

%kf
load('output_kf_run1.mat'); output_kf_run1 = output; clear output
load('output_kf_run2.mat'); output_kf_run2 = output; clear output
load('output_kf_run3.mat'); output_kf_run3 = output; clear output
load('output_kf_run4.mat'); output_kf_run4 = output; clear output
load('output_kf_run5.mat'); output_kf_run5 = output; clear output
load('output_kf_run6.mat'); output_kf_run6 = output; clear output


%% 2 make array together
% vf
output_vf.xres = [output_vf_run1.xres; output_vf_run2.xres; output_vf_run3.xres; output_vf_run4.xres; output_vf_run5.xres; output_vf_run6.xres];
output_vf.weightTest = [output_vf_run1.weightTest, output_vf_run2.weightTest, output_vf_run3.weightTest, output_vf_run4.weightTest, output_vf_run5.weightTest, output_vf_run6.weightTest]';
output_vf.errorData = [output_vf_run1.errorData; output_vf_run2.errorData; output_vf_run3.errorData; output_vf_run4.errorData; output_vf_run5.errorData; output_vf_run6.errorData];
output_vf.errorHaldane = [output_vf_run1.errorHaldane; output_vf_run2.errorHaldane; output_vf_run3.errorHaldane; output_vf_run4.errorHaldane; output_vf_run5.errorHaldane; output_vf_run6.errorHaldane];
output_vf.errorRegpars = [output_vf_run1.errorRegpars; output_vf_run2.errorRegpars; output_vf_run3.errorRegpars; output_vf_run4.errorRegpars; output_vf_run5.errorRegpars; output_vf_run6.errorRegpars];
output_vf.eData = [output_vf_run1.eData; output_vf_run2.eData; output_vf_run3.eData; output_vf_run4.eData; output_vf_run5.eData; output_vf_run6.eData];
output_vf.eHaldane = [output_vf_run1.eHaldane; output_vf_run2.eHaldane; output_vf_run3.eHaldane; output_vf_run4.eHaldane; output_vf_run5.eHaldane; output_vf_run6.eHaldane];
clear output_vf_run1 output_vf_run2 output_vf_run3 output_vf_run4 output_vf_run5 output_vf_run6

%kf
output_kf.xres = [output_kf_run1.xres; output_kf_run2.xres; output_kf_run3.xres; output_kf_run4.xres; output_kf_run5.xres; output_kf_run6.xres];
output_kf.weightTest = [output_kf_run1.weightTest, output_kf_run2.weightTest, output_kf_run3.weightTest, output_kf_run4.weightTest, output_kf_run5.weightTest, output_kf_run6.weightTest]';
output_kf.errorData = [output_kf_run1.errorData; output_kf_run2.errorData; output_kf_run3.errorData; output_kf_run4.errorData; output_kf_run5.errorData; output_kf_run6.errorData];
output_kf.errorHaldane = [output_kf_run1.errorHaldane; output_kf_run2.errorHaldane; output_kf_run3.errorHaldane; output_kf_run4.errorHaldane; output_kf_run5.errorHaldane; output_kf_run6.errorHaldane];
output_kf.errorRegpars = [output_kf_run1.errorRegpars; output_kf_run2.errorRegpars; output_kf_run3.errorRegpars; output_kf_run4.errorRegpars; output_kf_run5.errorRegpars; output_kf_run6.errorRegpars];
output_kf.eData = [output_kf_run1.eData; output_kf_run2.eData; output_kf_run3.eData; output_kf_run4.eData; output_kf_run5.eData; output_kf_run6.eData];
output_kf.eHaldane = [output_kf_run1.eHaldane; output_kf_run2.eHaldane; output_kf_run3.eHaldane; output_kf_run4.eHaldane; output_kf_run5.eHaldane; output_kf_run6.eHaldane];
clear output_kf_run1 output_kf_run2 output_kf_run3 output_kf_run4 output_kf_run5 output_kf_run6


%% 3 plot

%vf
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
sourceVm = 'experimentalSlopesFixed';
ode_pH = 'on';
setup.params = {'v_{maxFWD} [mM s^{-1}]'; 'K_{gap} [mM]'; 'K_{bpg} [mM]'; 'K_{nad} [mM]'; 'K_{nadh} [mM]'; 'v_{maxREV} [mM s^{-1}]'};
tempPars = output_vf.xres{16};
pvals = zeros(10,6);
plength = 6;
eQvals = [1.7E-3, 2.3E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1, 1.78E-1, 2.45E-1];
setup.pH_Keq_gapdh_eQ = 0.05 * (eQvals);
numpH = 10;
%     % Keq_pgk obtained with eQuilibrator
%     setup.pH_Keq_pgk = [1/(7.4E-4),    1/(7.2E-4),    1/(6.9E-4),    1/(6.5E-4),    1/(6.1E-4),    1/(5.7E-4),    1/(5.5E-4),    1/(5.3E-4),    1/(5.2E-4),   1/(5.2E-4)];
setup.exp_vmax_gapdhr = [0.001672560425949,   0.001708335558995,   0.001863340099700,   0.002134052317366,   0.001625999420026,   0.001392974025436,   0.001102411496380,   0.001004022996442,   0.000986761856102,   0.000877728986288];
setup.exp_vmax_gapdhf = [1.878012068989313e-04, 8.699614731347554e-05, 1.222088736070264e-04, 1.498266981509864e-04, 2.306088349420712e-04, 4.253144979769960e-04, 5.882596627863640e-04, 8.174876065012354e-04, 9.314111327450747e-04, 0.001119902785258];
pvals(:,1) = tempPars(1);
pvals(:,6) = tempPars(2);
pvals(:,2) = [tempPars(3); tempPars(7); tempPars(11); tempPars(15); tempPars(19); tempPars(23); tempPars(27); tempPars(31); tempPars(35); tempPars(39)];
pvals(:,3) = [tempPars(4); tempPars(8); tempPars(12); tempPars(16); tempPars(20); tempPars(24); tempPars(28); tempPars(32); tempPars(36); tempPars(40)];
pvals(:,4) = [tempPars(5); tempPars(9); tempPars(13); tempPars(17); tempPars(21); tempPars(25); tempPars(29); tempPars(33); tempPars(37); tempPars(41)];
pvals(:,5) = [tempPars(6); tempPars(10); tempPars(14); tempPars(18); tempPars(22); tempPars(26); tempPars(30); tempPars(34); tempPars(38); tempPars(42)];
pHvals = [6.1900    6.2600    6.4100    6.6000    6.8100    7.0600    7.2900    7.5100    7.6800    7.8100];
setup.pH_vals = [6.1900    6.2600    6.4100    6.6000    6.8100    7.0600    7.2900    7.5100    7.6800    7.8100];
% setup.exp_vmax_gapdhf(6); % mM s^{-1}
% setup.exp_vmax_gapdhr(6); % mM s^{-1}
errorData = [sum(abs(output_vf.errorData{16}(1:26))),...
            sum(abs(output_vf.errorData{16}(27:52))),...
            sum(abs(output_vf.errorData{16}(53:78))),...
            sum(abs(output_vf.errorData{16}(79:104))),...
            sum(abs(output_vf.errorData{16}(105:130))),...
            sum(abs(output_vf.errorData{16}(131:156))),...
            sum(abs(output_vf.errorData{16}(157:182))),...
            sum(abs(output_vf.errorData{16}(183:208))),...
            sum(abs(output_vf.errorData{16}(209:234))),...
            sum(abs(output_vf.errorData{16}(235:260)))];
errorDataRef = [sum(abs(output_vf.errorData{1}(1:26))),...
            sum(abs(output_vf.errorData{1}(27:52))),...
            sum(abs(output_vf.errorData{1}(53:78))),...
            sum(abs(output_vf.errorData{1}(79:104))),...
            sum(abs(output_vf.errorData{1}(105:130))),...
            sum(abs(output_vf.errorData{1}(131:156))),...
            sum(abs(output_vf.errorData{1}(157:182))),...
            sum(abs(output_vf.errorData{1}(183:208))),...
            sum(abs(output_vf.errorData{1}(209:234))),...
            sum(abs(output_vf.errorData{1}(235:260)))];

c = 'black';
figure(31)
for i = 1:plength
    
    % plot parameter values
    subplot(3,3,i)
    for o = 1
%         pvals = pvals_cell{o};
        plot(pHvals,pvals(:,i),'.-','color',c(o,:))
        hold on
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    end
    titleName = setup.params{i};
    title(titleName);
    
    % plot error Data
    if i == plength
        subplot(3,3,i+1)
        for o = 1
%             errorData = errorData_cell{o};
            plot(pHvals,errorData,'.-','color',c(o,:))
            hold on
            plot(pHvals,errorDataRef,'.-','color','blue')
        end
        title('error_{Data}')
        legend('error_{Data.wH1E3}','error_{Data.wH0}','Location','EastOutside')
    end
    
    % plot haldane relationship
    if i == plength
        Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
        Keq_haldane_estimated = zeros(1,numpH);
        for o = 1
%             pvals = pvals_cell{o};
            for j = 1:numpH
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
                switch ode_pH
                    case 'on'
                        H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
                    otherwise
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
            end
            subplot(3,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated,'.-','color',[.5 .5 .5])
            hold on
%             if o == 10
                semilogy(pHvals,Keq_haldane_theory,'ko','MarkerSize',5)
%             end
        end
%         legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
    end
    
end
suptitle('Km variable, Vm fixed')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  


%kf
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
sourceVm = 'experimentalSlopesFixed';
ode_pH = 'on';
setup.params = {'v_{maxFWD} [mM s^{-1}]'; 'K_{gap} [mM]'; 'K_{bpg} [mM]'; 'K_{nad} [mM]'; 'K_{nadh} [mM]'; 'v_{maxREV} [mM s^{-1}]'};
tempPars = output_kf.xres{12};
pvals = zeros(10,6);
plength = 6;
eQvals = [1.7E-3, 2.3E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1, 1.78E-1, 2.45E-1];
setup.pH_Keq_gapdh_eQ = 0.05 * (eQvals);
numpH = 10;
%     % Keq_pgk obtained with eQuilibrator
%     setup.pH_Keq_pgk = [1/(7.4E-4),    1/(7.2E-4),    1/(6.9E-4),    1/(6.5E-4),    1/(6.1E-4),    1/(5.7E-4),    1/(5.5E-4),    1/(5.3E-4),    1/(5.2E-4),   1/(5.2E-4)];
setup.exp_vmax_gapdhr = [0.001672560425949,   0.001708335558995,   0.001863340099700,   0.002134052317366,   0.001625999420026,   0.001392974025436,   0.001102411496380,   0.001004022996442,   0.000986761856102,   0.000877728986288];
setup.exp_vmax_gapdhf = [1.878012068989313e-04, 8.699614731347554e-05, 1.222088736070264e-04, 1.498266981509864e-04, 2.306088349420712e-04, 4.253144979769960e-04, 5.882596627863640e-04, 8.174876065012354e-04, 9.314111327450747e-04, 0.001119902785258];
pvals(:,2) = tempPars(1);
pvals(:,3) = tempPars(2);
pvals(:,4) = tempPars(3);
pvals(:,5) = tempPars(4);
pvals(:,1) = [tempPars(5); tempPars(7); tempPars(9); tempPars(11); tempPars(13); tempPars(15); tempPars(17); tempPars(19); tempPars(21); tempPars(23)];
pvals(:,6) = [tempPars(6); tempPars(8); tempPars(10); tempPars(12); tempPars(14); tempPars(16); tempPars(18); tempPars(20); tempPars(22); tempPars(24)];
pHvals = [6.1900    6.2600    6.4100    6.6000    6.8100    7.0600    7.2900    7.5100    7.6800    7.8100];
setup.pH_vals = [6.1900    6.2600    6.4100    6.6000    6.8100    7.0600    7.2900    7.5100    7.6800    7.8100];
% setup.exp_vmax_gapdhf(6); % mM s^{-1}
% setup.exp_vmax_gapdhr(6); % mM s^{-1}
errorData = [sum(abs(output_kf.errorData{12}(1:26))),...
            sum(abs(output_kf.errorData{12}(27:52))),...
            sum(abs(output_kf.errorData{12}(53:78))),...
            sum(abs(output_kf.errorData{12}(79:104))),...
            sum(abs(output_kf.errorData{12}(105:130))),...
            sum(abs(output_kf.errorData{12}(131:156))),...
            sum(abs(output_kf.errorData{12}(157:182))),...
            sum(abs(output_kf.errorData{12}(183:208))),...
            sum(abs(output_kf.errorData{12}(209:234))),...
            sum(abs(output_kf.errorData{12}(235:260)))];
errorDataRef = [sum(abs(output_kf.errorData{1}(1:26))),...
            sum(abs(output_kf.errorData{1}(27:52))),...
            sum(abs(output_kf.errorData{1}(53:78))),...
            sum(abs(output_kf.errorData{1}(79:104))),...
            sum(abs(output_kf.errorData{1}(105:130))),...
            sum(abs(output_kf.errorData{1}(131:156))),...
            sum(abs(output_kf.errorData{1}(157:182))),...
            sum(abs(output_kf.errorData{1}(183:208))),...
            sum(abs(output_kf.errorData{1}(209:234))),...
            sum(abs(output_kf.errorData{1}(235:260)))];

c = 'black';
figure(32)
for i = 1:plength
    
    % plot parameter values
    subplot(3,3,i)
    for o = 1
%         pvals = pvals_cell{o};
        plot(pHvals,pvals(:,i),'.-','color',c(o,:))
        hold on
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    end
    titleName = setup.params{i};
    title(titleName);
    
    % plot error Data
    if i == plength
        subplot(3,3,i+1)
        for o = 1
%             errorData = errorData_cell{o};
            plot(pHvals,errorData,'.-','color',c(o,:))
            hold on
            plot(pHvals,errorDataRef,'.-','color','blue')            
        end
        title('error_{Data}')
        legend('error_{Data.wH3E2}','error_{Data.wH0}','Location','EastOutside')
    end
    
    % plot haldane relationship
    if i == plength
        Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
        Keq_haldane_estimated = zeros(1,numpH);
        for o = 1
%             pvals = pvals_cell{o};
            for j = 1:numpH
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
%                     xassay(2) = x_temp(1);
%                     xassay(3) = x_temp(2);
%                     xassay(4) = x_temp(3);
%                     xassay(5) = x_temp(4);
%                         ks1 = 10 .^ x_temp(1) .* 2.48; % mM %k_gap
%                         ks2 = 10 .^ x_temp(3) .* 2.92; %mM %k_nad
%                         kp1 = 10 .^ x_temp(2) .* 1.18; % mM %k_bpg
%                         kp2 = 10 .^ x_temp(4) .* 0.022; % mM %k_nadh
% % % %                 Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                switch ode_pH
                    case 'on'
%                         H_effect = zeros(1:10);
%                         for j = 1:numpH
%                             H_effect(j) = 10^(setup.pH_vals(j) - setup.pH_vals(6));
%                         end
                        H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
    %                     Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                    otherwise
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
            end
            subplot(3,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated,'.-','color',[.5 .5 .5])
            hold on
%             if o == 10
                semilogy(pHvals,Keq_haldane_theory,'ko','MarkerSize',5)
%             end
        end
%         legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
    end
    
end
suptitle('Vm variable, Km fixed')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



%% 4 Save array
% 


%% 5 Save plot
% savefig(31,'fig03_parameter_overview_vf_wH1E3.fig');
% savefig(32,'fig03_parameter_overview_Kf_wH3E2.fig');

