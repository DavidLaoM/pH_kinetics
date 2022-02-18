% % Regularization + Results Visualization


%% Regularization
% lambdalist = [...
%     1E-5, 2E-5, 5E-5,...
%     1E-4, 2E-4, 5E-4,...
%     1E-3, 2E-3, 5E-3,...
%     1E-2, 2E-2, 5E-2,...
%     1E-1, 2E-1, 3E-1, 4E-1, 5E-1, 7E-1,... %area of change
%     1E0, 2E0, 3E0, 4E0, 5E0, 7E0,...%area of change
%     1E1, 2E1, 5E1,...
%     1E2, 2E2, 5E2,...
%     1E3, 2E3, 5E3,...
%     1E4, 2E4, 5E4,...
%     1E5, 2E5, 5E5,...
%     ];
% create array
eData = zeros(1,length(lambdalist));
eParameters = zeros(1,length(lambdalist));
for i = 1:length(lambdalist)
    eData(i) = sum(abs(array_eData{i}));
    eParameters(i) = sum(abs(array_eParams{i}));
end
% plotting
% f1 = figure(103);
yyaxis left
s1 = semilogx(lambdalist,eParameters,'o-','MarkerSize',6);
% ylabel('error_{Parameters} []')
% xlabel(erase(sprintf('selected lambda = %d',lambdalist(selLambdaPos)),".000000"))
hold on
yyaxis right
s2 = semilogx(lambdalist,eData,'o-','MarkerSize',6);
% ylabel('error_{Data} []')
l1 = line([lambdalist(selLambdaPos) lambdalist(selLambdaPos)],[s1.Parent.YLim(1) s1.Parent.YLim(2)]);
    l1.Color = 'black';
    l1.LineStyle = '--';
xlim([min(lambdalist) max(lambdalist)])
% set(103,'color','white')


