% here to be inserted in the main one
plPar_cell = cell(1,numpHtested);
plRes_cell = cell(1,numpHtested);
for j = 1:numpHtested
    fprintf('Start pH value #%d\n',j);
    x_temp = MSresult_DF12{i}.b;
    % select required data
    data.Vmaxs = data.Vmax(j,:);
    data.NADH = data.conc_mean(j,:);
    data.Vprofs = data.RRs(j,:);
    data.tempTime = data.time(j,:);

    limsPLA = setup.limsPLA;
    ptot = length(setup.params);
    names = setup.params;
    enzName = setup.enzymeName;
    
    for i = 1:ptot%(ptot-1)
        par         = x_temp; 
        index       = setup.caseStudy.parameters(i);
        parname     = setup.params{index};%name
        setup.i     = i;

        threshold   = chi2inv(0.5,size(par)); % threshold - chi square distribution, threshold = chi2inv(0.5,size(par));
        lb          = -limsPLA * ones(1,ptot); %lower bounds for parameter estimation (lsqnonlin)
        ub          = +limsPLA * ones(1,ptot); %upper bounds for parameter estimation (lsqnonlin)  

        func = @(x_temp)costfun_pH(x_temp,data,setup);
        Optimoptions= optimoptions('lsqnonlin','Display','off'); % options for optimization
        minStep     = 0.01; % minimal step factor
        maxStep     = 0.1; % maximal step factor
        minChange   = 0.001; % minimal change of resnorm
        maxChange   = 0.05; % maximal change of resnorm
        nr          = 100; % no. of samples in profile likelihood

        maxPar  = par(i) * 1.1;

        maxPar_up   = limsPLA;
        maxPar_down = -limsPLA;

        [plPar,plRes]=PLA_pH(func,par,i,maxPar,threshold,lb,ub,Optimoptions,minStep,maxStep,minChange,maxChange,nr,maxPar_up,maxPar_down,setup);
        pProfile.pPar.names{i}     = plPar';
        pProfile.pRes.names{i}     = plRes';
        fprintf('Finished: pH %f, parameter %f\n',j,i)
        saveName = sprintf('results/GAPDH_{reverse}_PLA_pH%d.mat',j);
        
        if i == ptot
            save(saveName,'pProfile');
        end

        if setup.plotResults == 1
            h = figure(j);

            subplot(3,3,i)
            plot(plPar, plRes)
            title([parname])
            xlabel('Parameter value')
            ylabel('Fit residual')

            xlim([-limsPLA limsPLA])

            hold on
            [val,ind] = min(plRes(~ismember(plRes,0)));
            ylim([min(plRes)-0.1 max(plRes)+0.1])
            y1 = ylim;
            line([par(i) par(i)],y1,'Color','magenta','LineStyle',':') % this is the right parameter value
            hold on
            line([plPar(ind) plPar(ind)],y1,'Color','red','LineStyle','--') % this is the current minimum
            hold on
%             drawnow()
        else
            h = 0;
        end
    end
    suptitleName = [enzName, ' PLA ', sprintf('pH %d', pHarray(j))];
    suptitle(suptitleName);
    
    plPar_cell{j} = plPar;
    plRes_cell{j} = plRes;
end