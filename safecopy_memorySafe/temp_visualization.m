% visualization
selLambdaPos = 1; % just adding a label to highlight in the map which lambda was chosen


% create array
eData = zeros(1,length(lambdalist));
eParameters = zeros(1,length(lambdalist));
for i = 1:length(lambdalist)
    eData(i) = sum(abs(array_eData{i}));
    eParameters(i) = sum(abs(array_eParams{i}));
end
% plotting
f1 = figure(103);
yyaxis left
s1 = semilogx(lambdalist,eParameters,'o-','MarkerSize',6);
ylabel('error_{Parameters} []')
xlabel(['regularization factor []: ',erase(sprintf('selected lambda = %d',lambdalist(selLambdaPos)),".000000")])
hold on
yyaxis right
s2 = semilogx(lambdalist,eData,'o-','MarkerSize',6);
ylabel('error_{Data} []')
textHere = {['Regulatization ',setup.enzymeName,': '];'Errors in parameter values and data fit in blue and red, ';'respectively, in the y-axes. Regularization factor in the x-axis'};
suptitle(textHere)
l1 = line([lambdalist(selLambdaPos) lambdalist(selLambdaPos)],[s1.Parent.YLim(1) s1.Parent.YLim(2)]);
    l1.Color = 'black';
    l1.LineStyle = '--';
set(103,'color','white')