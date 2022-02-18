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

%plot
figure(10)

%vf
subplot(1,2,1)
yyaxis left
semilogx(output_vf.weightTest,output_vf.eHaldane,'o-')
% loglog(output_vf.weightTest,output_vf.eHaldane,'o-')
ylabel('error_{Haldane.Keq}')
hold on
area([1E1 3E6],[35E-3 35E-3],'FaceColor',[0.9 0.9 0.9],'LineStyle','none');
set(gca,'children',flipud(get(gca,'children')))
hold on
yyaxis right
semilogx(output_vf.weightTest,output_vf.eData,'o-')
% loglog(output_vf.weightTest,output_vf.eData,'o-')
ylabel('error_{Data}')
legend('good fit region','e_{Haldane}', 'e_{Data}','Location','southoutside','orientation','horizontal')
title('Km variable, Vm fixed')
xlabel('weight_{Haldane}, weight_{Data} = 1')


%kf
subplot(1,2,2)
yyaxis left
semilogx(output_kf.weightTest, output_kf.eHaldane, 'o-')
% loglog(output_kf.weightTest, output_kf.eHaldane, 'o-')
ylabel('error_{Haldane.Keq}')
hold on
area([1E1 1E6],[1E-2 1E-2],'FaceColor',[0.9 0.9 0.9],'LineStyle','none');
set(gca,'children',flipud(get(gca,'children')))
hold on
yyaxis right
semilogx(output_kf.weightTest, output_kf.eData, 'o-')
% loglog(output_kf.weightTest, output_kf.eData, 'o-')
ylabel('error_{Data}')
legend('good fit region','e_{Haldane}', 'e_{Data}','Location','southoutside','orientation','horizontal')
title('Vm variable, Km fixed')
xlabel('weight_{Haldane}, weight_{Data} = 1')

suptitle('Haldane Regularization')


%% 4 Save array
% save('estimation_fixedParams.mat','output_vf','output_kf');


%% 5 Save plot
% savefig(10,'fig02_estimation_fixedParams.fig');

