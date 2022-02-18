function [PSAresult] = PSA(setup,data)
% function that carries out PSA over a specific enzyme. Used in the pH data
% study.
nlins = setup.nLinSpace;
legenda = setup.params;
linSpace = linspace(-1,1,nlins+1);
obsMet = setup.observableMetabolite;
ode = setup.ode;

if setup.enzymeName == 'gapdhr'
    setup.caseStudy = (1:7);
    switch ode
        case 'gapdhr_s_revMM'
            disp('Kinetics: simple reversible MM(vanEunen2012?)')
            disp('PHOS: not included');
        case 'vanHeerden2014'
            disp('Kinetics: as in van Heerden 2014 (GAPDH + PGK)')
            disp('PHOS: included');
        otherwise
            disp('WARNING: No ode has been selected');
    end
end
nparams = length(setup.caseStudy);

PSAsims = cell(nparams,nlins);
for i = 1:nparams
    for j = 1:nlins
        xtemp = zeros(1,nparams);
        xtemp(i) = linSpace(j);
%         data.chosenLink = data.DF(1,4);
        PSAsims{i,j} = simSys(xtemp,data,setup);
    end
end

if setup.plotResults == 1
    figure,
%     for i = 1:(nparams-1)
%         subplot(2,3,i)
    for i = 1:nparams
        subplot(3,3,i)
        for j = 1:nlins
            Tdata = PSAsims{i,j}.t;
            Ydata = PSAsims{i,j}.y(:,obsMet);
            plot(Tdata,Ydata,'color',[0.7 0.7 0.7])
            hold on
        end
        Tdata = PSAsims{i,setup.PSArefval}.t;
        Ydata = PSAsims{i,setup.PSArefval}.y(:,obsMet);
        plot(Tdata,Ydata,'color','b')
%         stdshade(data.raw.conc',0.1,'b',data.raw.time')
        title(legenda(i))
    end
    suptitleName = {'PSA GAPDH_{reverse}';'\color{blue}Literature parameter  \color[rgb]{0.7 0.7 0.7}Deviations */1'};
    suptitle(suptitleName);
end

[PSAresult] = PSAsims;
end