function [MSresult] = MS_pH(setup,data)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
pHtested = setup.pHtested;
numpHtested = nnz(pHtested);
limsMS = setup.limsMS;
nMS = setup.MSruns;
plength = length(setup.params);

x_temp = zeros(1,plength);
switch limsMS
    case 1
        ub = 1*ones(1,plength);
        ub(end) = 10;
        lb = -1*ones(1,plength);
        lb(end) = -10;
    otherwise
        disp('no boundaries specified in MultiStart');
end

% MS parameter estimation
MSresult = cell(numpHtested,1);
for i = 1:numpHtested    
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    model = @(x_temp)costfun_pH(x_temp,data,setup);
    problem = createOptimProblem('lsqnonlin', 'objective', model, ...
        'xdata', data, 'ydata', data, 'x0', x_temp, 'lb', lb, ...
        'ub', ub);
    ms = MultiStart('Display','iter');
    tic
    [b,fval,exitflag,output,solutions] = run(ms, problem, nMS);
    t = toc;
    fprintf('MS Pest finished for pH #%d, time %d s\n',i,t);
    res.b = b;
    res.fval = fval;
    res.exitflag = exitflag;
    res.output = output;
    res.solutions = solutions;
    res.t = t;
    MSresult{i} = res;
end

% plotting output histogram
if setup.plotResults == 1
    for k = 1:numpHtested 
%     for k = 1
        % recall the data
        [~,l2] = size(MSresult{k}.solutions);
        tempRes = zeros(l2,length(setup.params));
        for j = 1:l2
            tempRes(j,:) = MSresult{k}.solutions(j).X;
        end
        % plot data
        figure
        for j = 1:(length(setup.params))
            subplot(3,3,j)
            histogram(tempRes(:,j),setup.figWindow)
            title(setup.params{j})
%             set(gca,'xscale','log')
            xlim([-1 1])
        end
        pHarray = unique(data.pH);
        suptitleName = sprintf('MS histogram pH = %d',pHarray(k));
    %     newStr = erase(str,"the ")
        suptitle(suptitleName);
    end
end

end