% % fir09_plotKeqEffect
% In this plot the difference between the vm value estimated with the naïve
% and the computational approach is plotted. The biggest effect should be
% due to Keq change.
load('output_ald.mat');
load('output_eno.mat');
load('output_hxk.mat');
load('output_pgi.mat');
load('output_pyk.mat');

figure(1009)
%% subplot(2,3,1): GAPDH

%% subplot(2,3,2): ENO xres_selected(3:14)
subplot(2,3,2)
plot(output_eno.pHarray, output_eno.xres_selected(3:14),'k.-')
title('\color{blue}ENO')

%% subplot(2,3,3): HXK xres_selected(5:end)
subplot(2,3,3)
plot(output_hxk.pHarray, output_hxk.xres_selected(5:end),'k.-')
title('\color{red}HXK')

%% subplot(2,3,4): ALD xres_selected(4:11)
subplot(2,3,4)
plot(output_ald.pHarray, output_ald.xres_selected(4:11),'k.-')
title('\color{blue}ALD')

%% subplot(2,3,5): PYK xres_selected(7:18)
subplot(2,3,5)
plot(output_pyk.pHarray, output_pyk.xres_selected(7:18),'k.-')
title('\color{red}PYK')

%% subplot(2,3,6): PGI xres_selected(3:14)
subplot(2,3,6)
plot(output_pgi.pHarray, output_pgi.xres_selected(3:14),'k.-')
title('\color{blue}PGI')

%%
suptitle({'Change in Vm vs experimental approach (logarithmic)';...
    '\color{red} more \color{black} or \color{blue} less \color{black} change'})
%%
set(1009,'color','w');
savefig(1009,'fi09_KeqEffect');