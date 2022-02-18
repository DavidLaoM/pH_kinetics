% This follows the code and explanation Natal wrote in the pdf 
% 'nvr_Template for parameter estimation with Matlab Optimization Toolbox' 
% or in the following link 
% 'http://bmi.bmt.tue.nl/sysbio/parameter_estimation/Template%20for%20parameter%20estimation%20with%20Matlab%20Optimization%20Toolbox.pdf'.
% 
% All data required in:
% http://bmi.bmt.tue.nl/sysbio/parameter_estimation/parameter_estimation.html
%   
% First case

function compartment
%compartment Estimating pharmacokinetic parameters using a dynamic 2
%compartment model
%Loads datafile pmpat010.txt.
%% Load data
datafile = 'pmpat010.txt';
data = load(datafile);
t=data(1,:); %[s]
Ca=data(2,:); %[mM] arterial input function
Ct=data(3,:); %[mM] tissue response
figure; plot(t,Ca,'.b',t,Ct,'.r'); hold on
xlabel('Time [s]'); ylabel('[mM]')
legend('AIF','Ct')
N=length(Ct);
%% Simulate model with initial parameter values
%Simulate
Ktrans = 2e-3; ve = 0.1;
p = [Ktrans, ve];
[y,t] = compart_lsim(t,Ca,p);
plot(t,y,'k');
%% Estimate parameters and variances
optim_options = optimset('Display', 'iter',...
...%'TolFun', 1e-6,... %default: 1e-4
...%'TolX', 1e-6,... %default: 1e-4
'LevenbergMarquardt', 'on'); %default: 'off'
%optim_options = [];
p0 = [Ktrans, ve];
[p,resnorm,residual,exitflag,OUTPUT,LAMBDA,Jacobian] =...
lsqnonlin(@compart_error, p0, [],[],optim_options, t,Ca,Ct);
disp(' ')
p
%Estimated parameter variance
Jacobian = full(Jacobian); %lsqnonlin returns the Jacobian as a sparse
matrix
varp = resnorm*inv(Jacobian'*Jacobian)/N
stdp = sqrt(diag(varp)); %standard deviation is square root of variance
stdp = 100*stdp'./p; %[%]
disp([' Ktrans: ', num2str(p(1)), ' +/- ', num2str(stdp(1)), '%'])
disp([' ve: ', num2str(p(2)), ' +/- ', num2str(stdp(2)), '%']);
%% Simulate estimated model
[y,t] = compart_lsim(t,Ca,p);
figure; subplot(211); plot(t,Ct,'.r',t,y,'b');
xlabel('Time [s]'); ylabel('[mM]')
legend('data','model')