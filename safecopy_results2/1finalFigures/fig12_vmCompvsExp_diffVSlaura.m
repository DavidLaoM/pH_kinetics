% % fig12_vmCompvsExp_difffVSlaura
% In this file, the vm estimated via the computational method are compared
% to Laura's value, but in percentage.

% step 1. load
% step 2. recalculate
% step 3. plot
% (if needed) step 4. analyze


%% step 1. load
load('output_hxk.mat');
load('output_pgi.mat');
load('output_ald.mat');
% load('output_gapdhr.mat');
load('output_eno_kmfixed.mat');
load('output_pyk.mat');


%% step 2. recalculate
xdiffHXK = zeros(size(output_hxk)); % 12 dps
xselHXK = output_hxk.xres_selected(5:16);
xdiffPGI = zeros(size(output_pgi)); % 12 dps
xselPGI = output_pgi.xres_selected(3:14);
xdiffALD = zeros(size(output_ald)); % 8 dps
xselALD = output_ald.xres_selected(3:10);
xdiffENO = zeros(size(output_eno_kmfixed)); % 12 dps
xselENO = output_eno_kmfixed.xres_selected(3:14);
xdiffPYK = zeros(size(output_pyk)); % 12 dps
xselPYK = output_pyk.xres_selected(7:18);

for i = 1:12
    xdiffHXK(i) = ((10^xselHXK(i))-1)*100;
    xdiffPGI(i) = ((10^xselPGI(i))-1)*100;
    xdiffENO(i) = ((10^xselENO(i))-1)*100;
    xdiffPYK(i) = ((10^xselPYK(i))-1)*100;
end
for i = 1:8
    xdiffALD(i) = ((10^xselALD(i))-1)*100;
end


%% step 2.1. plot % difference
figure(34)

subplot(3,2,1) % vmHXK
plot(output_hxk.pHarray,xdiffHXK,'.-','color','black')
ylabel('v_{HXK}, %difference, NADPH')
ylim([-100 100])
line([6 8],[-25 -25],'color','red','linestyle','--')
line([6 8],[25 25],'color','red','linestyle','--')
line([6 8],[0 0],'color','blue','linestyle','--')

subplot(3,2,2) % vmPGI
plot(output_pgi.pHarray,xdiffPGI,'.-','color','black')
ylabel('v_{PGI}, %difference, NADH')
ylim([-100 100])
line([6 8],[-25 -25],'color','red','linestyle','--')
line([6 8],[25 25],'color','red','linestyle','--')
line([6 8],[0 0],'color','blue','linestyle','--')

subplot(3,2,3) % vmALD
plot(output_ald.pHarray,xdiffALD,'.-','color','black')
ylabel('v_{ALD}, %difference, NADH')
ylim([-100 100])
line([6 8],[-25 -25],'color','red','linestyle','--')
line([6 8],[25 25],'color','red','linestyle','--')
line([6 8],[0 0],'color','blue','linestyle','--')

% vmGAPDHR

subplot(3,2,5) % vmENO
plot(output_eno_kmfixed.pHarray,xdiffENO,'.-','color','black')
ylabel('v_{ENO}, %difference, PEP')
ylim([-100 100])
line([6 8],[-25 -25],'color','red','linestyle','--')
line([6 8],[25 25],'color','red','linestyle','--')
line([6 8],[0 0],'color','blue','linestyle','--')

subplot(3,2,6) % vmPYK
plot(output_pyk.pHarray,xdiffPYK,'.-','color','black')
ylabel('v_{PYK}, %difference, NADH')
% ylim([-100 100])
line([6 8],[-25 -25],'color','red','linestyle','--')
line([6 8],[25 25],'color','red','linestyle','--')
line([6 8],[0 0],'color','blue','linestyle','--')

suptitleText = {'v_{m} [umol_{metaboliteAssay} mg_{P}^{-1} min^{-1}]';'Not normalized';'\color{blue}Laura.s values \color{black}and \color{red}25% deviation limits'};
suptitle(suptitleText);
set(gcf,'color','w');


