% % load the data (if needed. Better not, since it will change when
% % overwritting)
% loadName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_regularizationResults.mat'];
% load(loadName);
% create blank arrays
n = length(lambdalist);
cell_N = cell(1,n);
cell_Jacobian = cell(1,n);
cell_varp = cell(1,n);
cell_stdp = cell(1,n);
% fill in the arrays
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
for i = 1:n
    % recall
    setup.selectedLambda = lambdalist(i);
    tempError = array_raw_error{i};
    tempJacobian = array_Jacobian{i};
    resnorm = array_resnorm{i};
    % calculate
    N = length(tempError);
    Jacobian = full(tempJacobian);  
    varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
    stdp = sqrt(diag(varp));
    % save
    cell_N{i} = N;
    cell_Jacobian{i} = Jacobian;
    cell_varp{i} = varp;
    cell_stdp{i} = stdp;
end

% pre treatment for the plotting
m = length(lambdalist);
n = length(setup.parameterNames);
parameterEstimate = zeros(n,m);
confidenceInterval = zeros(n,m);
for i = 1:n
    for j = 1:m
        parameterEstimate(i,j) = array_xres{j}(i);
        confidenceInterval(i,j) = cell_stdp{j}(i);
    end
end

% plotting
figure(104)
for i = 1:n
    sp = subplot(5,5,i);
    yyaxis left
    semilogx(lambdalist,parameterEstimate(i,:),'.-')
    if((i == 1)||(i == 6)||(i == 11)||(i == 16)||(i == 21))
%         ylabel({'Parameter Estimate:'; setup.parameterNames{i}},'FontSize',10)
        ylabel('Parameter Estimate','FontSize',10)
    end
    if(n<=5)
        n_idxs = [1 2 3 4 5];
    elseif((n>=6)&&(n<=10))
        n_idxs = [6 7 8 9 10];
    elseif((n>=11)&&(n<=15))
        n_idxs = [11 12 13 14 15];
    elseif((n>=16)&&(n<=20))
        n_idxs = [16 17 18 19 20];
    elseif(n>=21)
        n_idxs = [21 22 23 24 25];
    end
    if((i == n_idxs(1))||(i == n_idxs(2))||(i == n_idxs(3))||(i == n_idxs(4))||(i == n_idxs(5)))
        xlabel('regularization factor','FontSize',10)
    end
    yyaxis right
    semilogx(lambdalist,confidenceInterval(i,:),'.-')
    if((i == 5)||(i == 10)||(i == 15)||(i == 20)||(i == 25))
%         ylabel({'Confidence interval:'; setup.parameterNames{i}},'FontSize',10)
        ylabel('Confidence interval:','FontSize',10)
    end
    title(setup.parameterNames{i},'FontSize',10)
    xlim([min(lambdalist) max(lambdalist)])
    % lambdalistbar
    l1 = line([lambdalist(selLambdaPos) lambdalist(selLambdaPos)],[sp.YLim(1) sp.YLim(2)]);
    l1.Color = 'black';
    l1.LineStyle = '--';
end
textHere = [setup.enzymeName,': parameters and confidence intervals vs regularization factor'];
suptitle(textHere)
set(104,'color','white')

