% % expRatesDetermination
% 

% blank cell where the rates are stored, also taken time, start point and
% r2 coefficient
Vmax_mw = cell(size(dp_start));
Vmax_mw_start = cell(size(dp_start));
num_mw = cell(size(dp_start));
Vmax_mw_time = cell(size(dp_start));
R2_mw = cell(size(dp_start));
% cell where the vmax is selected
idxs_opt_pre = zeros(size(dp_start));
idxs_opt = zeros(size(dp_start));
Vmax_mw_opt = zeros(size(dp_start));
Vmax_mw_opt_corr = zeros(size(dp_start));
R2_mw_opt = zeros(size(dp_start));

% loop for calculating Vmax_mw and Vmax_mw_opt
for i = 1:numpHtested
    for j = 1:DFs
        if setup.caseStudyGAPDH == 1
            for case1 = 1
                % method 1
                if(j <= 3)
                    minwindow = 40;
                elseif((j >= 4)&&(j <= 7)&&(j == 12)) 
                    minwindow = 20;
                elseif((j >= 8)&&(j <= 11))
                    minwindow = 10;
                end
                % method 2
                if(j >= 4)
                    minwindow = 5;
                end
                % method 3
                if(i >= 2)
                    minwindow = 40;
                end
                % method 3 (2020 - 10 - 22)
                minwindowVals = ones(1,numpHtested);
                minwindowVals(1) = 30; %
                minwindowVals(2) = 5;
                minwindowVals(3) = 30;
                minwindowVals(4) = 8;
                minwindowVals(5) = 6; %
                minwindowVals(6) = 4;
                minwindowVals(7) = 6;
                minwindowVals(8) = 4;
                minwindowVals(9) = 5; %
                minwindowVals(10) = 5;
                minwindowVals(11) = 4;
                minwindowVals(12) = 8;
                minwindow = minwindowVals(i);
            end
        elseif setup.caseStudyGAPDHr == 1
            for case1 = 1
                if j == 4
                    minwindow = 3;
                else
                    minwindow = 6;
                end
            end
        elseif setup.enzymeName == 'pdc'
            for case1 = 1
                if((j == 4)&&(i <= 6))
                    minwindow = 3;
                elseif j == 3
                    minwindow = 15;
                elseif j == 2
                    minwindow = 40;
                elseif j == 1
                    minwindow = 40;
                else
                    minwindow = 6;
                end
            end
        elseif setup.enzymeName == 'pgm'
            for case1 = 1
                if((i == 2)||(i == 10))
                    minwindow = 4;
                else
                    minwindow = 7;
                end
            end
        elseif setup.enzymeName == 'pyk'
            for case1 = 1
                if((j == 4)&&(i >= 3))
                    minwindow = 4; 
                elseif((j == 4)&&(i <= 2))
                    minwindow = 1; 
                elseif((j == 3)&&(i <= 2))
                    minwindow = 10; 
                elseif((j == 3)&&(i >= 3))
                    minwindow = 20; 
                elseif((j == 2)&&(i <= 2))
                    minwindow = 30; 
                elseif((j == 2)&&(i >= 3))
                    minwindow = 50; 
                elseif((j == 1)&&(i <= 2))
                    minwindow = 50; 
                elseif((j == 1)&&(i >= 3))
                    minwindow = 75; 
                end
            end
        else
        end
        total_len(i,j) = length(data.conc_mean{i,j});
        Vmax_mw_temp = [];
        Vmax_mw_start_temp = [];
        R2_mw_temp = [];
        endVal = total_len(i,j)-1;
        for k = dp_start(i,j):endVal
            k2 = k+1;
            % calculate Vmax_mw
            tempVmax = zeros((total_len(i,j) - dp_start(i,j)),1);
            tempX = data.time{i,j}(dp_start(i,j):k2);
            tempY = data.conc_mean{i,j}(dp_start(i,j):k2);
            p = polyfit(tempX,tempY,1);
            if setup.caseStudyGAPDH == 1 %gapdh_fwd
                Vmax_mw_temp = [Vmax_mw_temp;p(1)];
            elseif setup.caseStudyGAPDHr == 1 %gapdh_fwd 
                Vmax_mw_temp = [Vmax_mw_temp;-p(1)];               
            elseif setup.enzymeName == 'hxk'
                Vmax_mw_temp = [Vmax_mw_temp;p(1)];
            elseif setup.enzymeName == 'eno'
                Vmax_mw_temp = [Vmax_mw_temp;p(1)];
            else
                Vmax_mw_temp = [Vmax_mw_temp;-p(1)];
            end
            % calculate Vmax_mw_temp
            Vmax_mw_start_temp = [Vmax_mw_start_temp;p(2)];
            % calculate R2
            yfit = polyval(p,tempX);
            yresid = tempY - yfit;
            SSresid = sum(yresid.^2);
            SStotal = (length(tempY)-1) * var(tempY);
            rsq = 1 - SSresid/SStotal;
            R2_mw_temp = [R2_mw_temp; rsq];
        end
        Vmax_mw_time{i,j} = data.time{i,j}(dp_start(i,j):endVal);
        Vmax_mw{i,j} = Vmax_mw_temp;
        Vmax_mw_start{i,j} = Vmax_mw_start_temp;
        R2_mw{i,j} = R2_mw_temp;
        num_mw_temp = 1:1:length(Vmax_mw_time{i,j});
        num_mw{i,j} = num_mw_temp';
        
        % select idxs optimum based on maximum, after minwindow points
        [~,idxs_opt_pre(i,j)] = max(Vmax_mw{i,j}(minwindow:end));
        idxs_opt(i,j) = idxs_opt_pre(i,j) + minwindow - 1;
        % specific for hxk
        if((setup.caseStudyGAPDH == 1)||(setup.caseStudyGAPDHr == 1)) %gapdh_fwd
        elseif setup.enzymeName == 'hxk'
            idxs_opt(i,j) = idxs_one_initial(i,j);
        end
        % calculate Vmax_mw_opt
        Vmax_mw_opt(i,j) = Vmax_mw{i,j}(idxs_opt(i,j));
        Vmax_mw_opt_corr(i,j) = Vmax_mw_opt(i,j) .* DF(i,j);
        R2_mw_opt(i,j) = R2_mw{i,j}(idxs_opt(i,j));
        
        % display exit
%         formatSpec = 'i = %f, j = %f\n';
%         fprintf(formatSpec,i,j)
    end
end

% %%
% dbstop if error
% select the dilution rate to consider
considerExpDF = zeros(size(dp_start));
for i = 1:numpHtested
    tempMean = mean(Vmax_mw_opt_corr(i,:));
    tempStd = std(Vmax_mw_opt_corr(i,:));
    tempRage = [tempMean-tempStd tempMean+tempStd];
    for j = 1:DFs
        tempVal = Vmax_mw_opt_corr(i,j);
        if((tempVal >= tempRage(1))&&(tempVal <= tempRage(2)))
            considerExpDF(i,j) = 1;
        else
        end
    end
end

temp_Vmax_mw_opt = Vmax_mw_opt; % safecopy for plot
temp_Vmax_mw_opt_corr = Vmax_mw_opt_corr; % safecopy for plot

% %% (0.2) plot moving window
% converting units
conversionFactor =  60 .* 60 ./ setup.concProtein;
Vmax_mw_opt = temp_Vmax_mw_opt * conversionFactor;
Vmax_mw_opt_corr = temp_Vmax_mw_opt_corr * conversionFactor * setup.branchFactor;

% Final determination of the experimental vmax
Vmax_experimental = zeros(size(pHarray));
stDev_experimental = zeros(size(pHarray));
for i = 1:length(Vmax_experimental)
    temp_vmax1 = Vmax_mw_opt_corr(i,:);
    temp_idx1 = idxs2consider(i,:);
    temp_vmax2 = [];
    for j = 1: length(temp_vmax1)
        if temp_idx1(j) == 1
            temp_vmax2 = [temp_vmax2; temp_vmax1(j)];
        end
    end
    Vmax_experimental(i) = mean(temp_vmax2);
    stDev_experimental(i) = std(temp_vmax2);
end

% figures
if setup.plotOutput == 1
    % vmaxs vs DF
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:numpHtested
        h1 = subplot(3,4,i);
        % values at different DFs
        yyaxis left, plot(DFarray,Vmax_mw_opt(i,:),'o-')
        h1.YLim = [0 h1.YLim(2)];
        % ylabel left
        if((i == 1)||(i == 5)||(i == 9))
            ylabel({'Vmax_{uncorrected}';'[umol mg_{P}^{-1} min^{-1}]'})
        end
        % xlabel
        if((setup.caseStudyALD == 1)||(setup.caseStudyTPI == 1))
            i_xlabel = 5:8;
        else
            i_xlabel = 9:12;
        end
        if((i == i_xlabel(1))||(i == i_xlabel(2))||(i == i_xlabel(3))||(i == i_xlabel(4)))
            xlabel('Dilution factor []')
        end
        hold on
        % diagonal line
        testReg = idxs2consider(i,:);
        o = 5;
        while o > 1
            o = o - 1;
            if testReg(o) == 1
                break
            end
        end
        line([0 DFarray(o)],[0 Vmax_mw_opt(i,o)],'Color','black','LineStyle','--')
        % corrected values
        hold on
        yyaxis right, plot(DFarray,Vmax_mw_opt_corr(i,:),'o-')
        for j = 1:DFs
            if any(idxs2consider(i,j)) == 1
                plot(DFarray(j), Vmax_mw_opt_corr(i,j)+2E-4, 'ko-','MarkerFaceColor','k')
            end
        end
        ylim([0.8*min(Vmax_mw_opt_corr(i,:)) 1.2*max(Vmax_mw_opt_corr(i,:))])
        % ylabel right
        if((i == 4)||(i == 8)||(i == 12))
            ylabel({'Vmax_{corrected}';'[umol mg_{P}^{-1} min^{-1}]'})
        end
        % horizontal line    
        hold on
        line(h1.XLim,[Vmax_experimental(i) Vmax_experimental(i)],'Color','black','LineStyle','--')   
        % title
        title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
    end
    suptitleName = [setup.enzymeName, ': Vmax vs DF after moving window, raw (blue) and corrected (red)'];
    suptitle(suptitleName)

    % vmax moving window
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:numpHtested
        subplot(3,4,i)
        for j = 1:DFs
            line([idxs_opt(i,j) idxs_opt(i,j)],limRates,'Color','black','LineStyle',':')
            hold on
            plot(num_mw{i,j},Vmax_mw{i,j},'.-')
            hold on
        end
        title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
    %     if i == numpHtested
    %         if setup.caseStudyGAPDHr == 1
    %             legend('DF 8','DF 4','DF 2','DF 1')
    %         end
    %     end
        % lim y-axis
        ylim(limRates)
    end
    suptitleName = ['Vmax moving window: ',setup.enzymeName];
    suptitle(suptitleName)

    % R2 moving window
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:numpHtested
        subplot(3,4,i)
        for j = 1:DFs
            line([idxs_opt(i,j) idxs_opt(i,j)],limR2,'Color','black','LineStyle',':')
            hold on
            plot(num_mw{i,j},R2_mw{i,j},'k.-')
            hold on
            text(num_mw{i,j}(end)+5,R2_mw{i,j}(end),sprintf('DF%d',DF(i,j)))
            hold on
        end
        % ylabel left
        if((i == 1)||(i == 5)||(i == 9))
            ylabel('R2 []')
        end
        % xlabel
        if((setup.caseStudyALD == 1)||(setup.caseStudyTPI == 1))
            i_xlabel = 5:8;
        else
            i_xlabel = 9:12;
        end
        if((i == i_xlabel(1))||(i == i_xlabel(2))||(i == i_xlabel(3))||(i == i_xlabel(4)))
            xlabel('moving window size []')
        end
        % title
        title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
    end
    suptitleName = ['R2 moving window: ',setup.enzymeName];
    suptitle(suptitleName)

    % vmax start moving window
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:numpHtested
        subplot(3,4,i)
        for j = 1:DFs
            line([idxs_opt(i,j) idxs_opt(i,j)],limcConc,'Color','black','LineStyle',':')
            hold on
            plot(num_mw{i,j},Vmax_mw_start{i,j},'.-')
            hold on
        end
        title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
    end
    suptitleName = ['Starting concentration moving window: ',setup.enzymeName];
    suptitle(suptitleName)

    % experimental vmax vs pH
    figure
    errorbar(pHarray,Vmax_experimental,stDev_experimental,'o-')
    suptitleName = ['Experimental Vmax vs pH: ',setup.enzymeName];
    suptitle(suptitleName)
end

% adjusting for the experimentally determined vmax (better before in the
% pipeline)
data.Vmax(:,4) = Vmax_experimental / (60 * 60 / setup.concProtein);

