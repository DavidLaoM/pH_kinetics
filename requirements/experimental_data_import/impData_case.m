function [output] = impData_case(setup)
% % gapdhR_impData.m
% import the dataset from the excel files from the gapdh_reverse reaction.

% steps
% 1- Determine the limits in the dataset sheet to import.
% 2- Select the specific ranges where data is found.
% 3- Import data. Organise in cells. Different layers for property/feature.
% 4- Calculation of mean, std and reaction rates
% 5- Visualization 
% 6- Saving output variables and figures

tempText = ['Starting import for enzyme ',setup.enzymeName];
f = waitbar(0,tempText);
pause(.5)

% % % initial pop up messages with script scope
% % disp('For this script to work:');
% % disp('Units should be kept the same as well');
% % disp('Excel sheets have to keep the same structure as gapdh_r.');
% % disp('All data points in one sheet should have the same time points.');
% % disp('In the plate reader, all rows in one column have to be associated to the same excel sheet.');
% % disp('WARNING: baseline experiments did not finish at the same time than the experiment itself.');
% % disp('WARNING: hxk datasets had the lowest sheet name replaced, from sheet 3 and 4 to 7 and 8, respectively. Addressed in the layout excel sheet as well')
% % disp('WARNING: pgi repetition, plate_3, is missed')

% (1) Determine the limits in the dataset sheet to import.
tempText = [setup.enzymeName, ' Step 1 of 6: Determining dataset locations'];
waitbar(1/6,f,tempText);
pause(1)
layoutFileName = setup.layoutFileName;
sheetNumberMax = setup.layoutSheetNumberMax;
sizePlates = setup.sizePlates;

tempNum = flip(1:sheetNumberMax);
for i = tempNum
    tempSheet = sprintf('plate_%d',i);
    [~,tempB,~] = xlsfinfo(layoutFileName);
    sheetValid = any(strcmp(tempB,tempSheet));
    if sheetValid == 1
       sheetNumber = i;
       break
    end
end

    % TO CLEAN UP missing PGI repetition
    if setup.caseStudyPGI == 1
        sheetNumber = 2;
    end

Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
corrAlph = Alphabet(strfind(Alphabet,setup.sizePlates(1)):end);
limits = cell(sheetNumber,4);
sizeData = zeros(sheetNumber,2);
for i = 1:sheetNumber
    sheetID = sprintf('plate_%d',i);
    [num,~,raw] = xlsread(layoutFileName,sheetID,sizePlates);
    sizeNum = size(num);
    raw2 = cell2mat(raw);
    raw3 = ~isnan(raw2);
    sizeRaw3 = size(raw3);
    [row,col] = find(raw3 == 1);
    
    limits{i,1} = corrAlph(col(1)); %letterStart
    limits{i,2} = corrAlph(col(1) + sizeNum(2) - 1); %letterEnd
    limits{i,3} = 2 + row(1) - 1; %numStart
    limits{i,4} = 2 + row(1) - 1 + sizeNum(1) - 1; %numEnd    
end
% disp(limits);

% (2) Select the specific ranges where data is found.
tempText = [setup.enzymeName, ' Step 2 of 6: Selecting dataset ranges'];
waitbar(2/6,f,tempText);
pause(1)
% feature = {'excelSheet' 'coordinates' 'pH' 'dilution'};
importRanges = cell(sheetNumber+1,4);
importRanges{end,1} = 'excelSheet';
importRanges{end,2} = 'coordinates';
importRanges{end,3} = 'pH';
importRanges{end,4} = 'dilution';
for i = 1:sheetNumber
    importRanges{i,1} = [limits{i,1},sprintf('%d:',limits{i,3}),limits{i,2},sprintf('%d',limits{i,4})];
    importRanges{i,2} = [limits{i,1},sprintf('%d:',limits{i,3}+10),limits{i,2},sprintf('%d',limits{i,4}+10)];
    importRanges{i,3} = [limits{i,1},sprintf('%d:',limits{i,3}+20),limits{i,2},sprintf('%d',limits{i,4}+20)];
    importRanges{i,4} = [limits{i,1},sprintf('%d:',limits{i,3}+30),limits{i,2},sprintf('%d',limits{i,4}+30)];
end
% % disp('Layout data selected from the following ranges')
% % disp(importRanges);

% (3) Import data. Organise in cells. Different layers for property/feature.
tempText = [setup.enzymeName, ' Step 3 of 6: Importing data'];
waitbar(3/6,f,tempText);
pause(1)
%layout data
%layout data: separate
tempExcelSheet = cell(sheetNumber,1);
tempCoordinates = cell(sheetNumber,1);
tempPH = cell(sheetNumber,1);
tempDilution = cell(sheetNumber,1);
for i = 1:sheetNumber
    sheetID = sprintf('plate_%d',i);
    [tempExcelSheet{i,1},~,~] = xlsread(layoutFileName,sheetID,importRanges{i,1});
    [~,tempCoordinates{i,1},~] = xlsread(layoutFileName,sheetID,importRanges{i,2});
%     [tempPH{i,1},~,~] = xlsread(layoutFileName,sheetID,importRanges{i,3});
    % changed since some pH values in the excel file are not annotated as
    % numeric, but string.
    [~,~,tempRaw] = xlsread(layoutFileName,sheetID,importRanges{i,3});
    [s1,s2] = size(tempRaw);
    tempTemp = cell(s1,s2);
    for j = 1:(s1*s2)
        stringORnum = isnumeric(tempRaw{j});
        if stringORnum == 1
            tempTemp{j} = tempRaw{j};
        elseif stringORnum == 0
            tempTemp{j} = str2double(tempRaw{j});
        end
    end    
    tempPH{i,1} = tempTemp;
    
    [tempDilution{i,1},~,~] = xlsread(layoutFileName,sheetID,importRanges{i,4});
end
% disp(tempExcelSheet);
% layout data: put together
excelSheet = [];
coordinates = [];
pH = [];
Dilution = [];
for i = 1:sheetNumber
    excelSheet = [excelSheet, tempExcelSheet{i}];
    coordinates = [coordinates, tempCoordinates{i}];
    pH = [pH, tempPH{i}];
    Dilution = [Dilution, tempDilution{i}];
end
% disp(excelSheet);

% progression curves
FileName = setup.FileName;
dataSheetNumberMax = setup.dataSheetNumberMax;
dataPointsMax = setup.dataPointsMax;

tempNum = 2*flip(1:dataSheetNumberMax);
for i = tempNum
    tempSheet = sprintf('Sheet%d',i);
    [~,tempB,~] = xlsfinfo(FileName);
    sheetValid = any(strcmp(tempB,tempSheet));
    if sheetValid == 1
       dataSheetNumber = i/2;
       break
    end
end

    % TO CLEAN UP missing PGI repetition
    if setup.caseStudyPGI == 1
        dataSheetNumber = 6;
    end


tempNum2 = 2*(1:dataSheetNumberMax);
for i = tempNum2
    tempSheet = sprintf('Sheet%d',i);
    [~,tempB,~] = xlsfinfo(FileName);
    sheetValid = any(strcmp(tempB,tempSheet));
    if sheetValid == 1
       dataSheetNumberFirst = i/2;
       break
    end
end

dataPoints = dataPointsMax;

lmax = setup.lettermax;
fileNameB = setup.fileNameBackground;

rawData = cell(dataSheetNumber - dataSheetNumberFirst + 1,1);
for i = dataSheetNumberFirst:dataSheetNumber
    sheetID = sprintf('Sheet%d',i*2);
%     disp(sheetID);
%     [~,~,rawData{i}] = xlsread(FileName,sheetID);
    [~,rawData{i},~] = xlsread(FileName,sheetID);
end
limits2 = size(excelSheet);
Alphabet2 = cell(lmax,length(Alphabet));
for i = 1:lmax
    for j = 1:length(Alphabet)
        if i == 1
            Alphabet2{i,j} = [Alphabet(j)];
        else
            Alphabet2{i,j} = [Alphabet(i-1),Alphabet(j)];
        end
    end
end
% tempSize = size(Alphabet2);
% Alphabet3 = cell(1,tempSize(1)*tempSize(2));
Alphabet3 = cell(1,1);
for i = 1:lmax
    Alphabet3 = {Alphabet3{1,:}, Alphabet2{i,:}};
    if i == 1
        Alphabet3(:,1) = [];
    end
end
% Alphabet3 = {Alphabet2{1,:}, Alphabet2{2,:}, Alphabet2{3,:}};

% absorbance_raw
absorbance_raw = cell(limits2);
for i = 1:limits2(1)
    for j =1:limits2(2)
        sheetID = sprintf('Sheet%d',excelSheet(i,j));
        tempID = strfind(rawData{excelSheet(i,j)/2},coordinates{i,j});
        temp1 = find(~cellfun(@isempty,tempID));
        temp2 = temp1(1);
        label = Alphabet3{temp2};
        rangeID = [label,sprintf('%d:',2),label,sprintf('%d',1+dataPoints)];
        absorbance_raw{i,j} = xlsread(FileName,sheetID,rangeID);
%         disp(i); disp(j);
    end
end

% % % WARNING
% % disp('Warning: the import of background data works as long as the layout between the experimenta and baseline file are the same.');

% absorbance_background
rawDataB = cell(dataSheetNumber - dataSheetNumberFirst + 1,1);
for i = dataSheetNumberFirst:dataSheetNumber
    sheetID = sprintf('Sheet%d',i*2-1);
    [~,rawDataB{i},~] = xlsread(fileNameB,sheetID);
end
absorbance_background = cell(limits2);
excelSheet2 = excelSheet - 1;
for i = 1:limits2(1)
    for j =1:limits2(2)
        sheetID = sprintf('Sheet%d',excelSheet2(i,j));
        tempID = strfind(rawDataB{(excelSheet2(i,j)+1)/2},coordinates{i,j});
        temp1 = find(~cellfun(@isempty,tempID));
        temp2 = temp1(1);
        label = Alphabet3{temp2};
        rangeID = [label,sprintf('%d:',2),label,sprintf('%d',1+dataPoints)];
        absorbance_background{i,j} = xlsread(fileNameB,sheetID,rangeID);
%         disp(i); disp(j);
    end
end

% time
time = cell(dataSheetNumber - dataSheetNumberFirst + 1,1);
time_background = cell(dataSheetNumber - dataSheetNumberFirst + 1,1);
for i = dataSheetNumberFirst:dataSheetNumber
    sheetIDe = sprintf('Sheet%d',i*2);
    [tempTime,~,~] = xlsread(FileName,sheetIDe,sprintf('A2:A%d',dataPoints)); % in seconds
    time{i} = 24*60*tempTime;
    sheetIDb = sprintf('Sheet%d',i*2-1);
    [tempTime_background,~,~] = xlsread(fileNameB,sheetIDb,sprintf('A2:A%d',dataPoints)); % in seconds
    time_background{i} = 24*60*tempTime_background;
end

% (4) Calculation of mean, std, concentration and reaction rates
tempText = [setup.enzymeName, ' Step 4 of 6: Calculations'];
waitbar(4/6,f,tempText);
pause(1)

extinction_coefficient = setup.extinction_coefficient;
nreps = setup.numberReplicates;
cic = setup.casesInColums;
delWrong = setup.deleteWrongDatapoints;

% absorbance_corrected
absorbance_background_interp = cell(limits2);
background = cell(limits2);
absorbance_corrected = cell(limits2);
for i = 1:limits2(1)
    for j =1:limits2(2)
        % added in case that the background array stopped before than the
        % experiment
        
        if((time_background{(excelSheet2(i,j)+1)/2}(end)+1) <= time{excelSheet(i,j)/2}(end))
            temp_time_background = [time_background{(excelSheet2(i,j)+1)/2}; time{excelSheet(i,j)/2}(end)];
            addedValue = absorbance_background{i,j}(end) + (temp_time_background(end) - temp_time_background(end-1))*((absorbance_background{i,j}(end)-absorbance_background{i,j}(1))/(temp_time_background(end-1)-temp_time_background(1)));
            temp_absorbance_background = [absorbance_background{i,j}; addedValue];
            absorbance_background_interp{i,j} = interp1(temp_time_background,temp_absorbance_background,time{excelSheet(i,j)/2},'pchip');
        else
            absorbance_background_interp{i,j} = interp1(time_background{(excelSheet2(i,j)+1)/2},absorbance_background{i,j},time{excelSheet(i,j)/2},'pchip');
        end
% % % %         absorbance_background_interp{i,j} = interp1(temp_time_background,temp_absorbance_background,time{excelSheet(i,j)/2},'pchip');
%         if(time_background{(excelSheet2(i,j)+1)/2}(end) <= time{excelSheet(i,j)/2}(end))
%             time_background{(excelSheet2(i,j)+1)/2} = [time_background{(excelSheet2(i,j)+1)/2}; time{excelSheet(i,j)/2}(end)];
%             addedValue = absorbance_background{i,j}(end) + (time_background{(excelSheet2(i,j)+1)/2}(end) - time_background{(excelSheet2(i,j)+1)/2}(end-1))*((absorbance_background{i,j}(end)-absorbance_background{i,j}(1))/(time_background{(excelSheet2(i,j)+1)/2}(end-1)-time_background{(excelSheet2(i,j)+1)/2}(1)));
%             absorbance_background{i,j} = [absorbance_background{i,j}; addedValue];
%         end 
%         absorbance_background_interp{i,j} = interp1(time_background{(excelSheet2(i,j)+1)/2},absorbance_background{i,j},time{excelSheet(i,j)/2},'pchip');
        background{i,j} = absorbance_background_interp{i,j} - absorbance_background_interp{i,j}(end);
        absorbance_corrected{i,j} = absorbance_raw{i,j} - background{i,j};
%         disp(i);disp(j);
    end
end

% absorbance mean and standard deviation
absorbance_mean = cell(limits2(1)/nreps,limits2(2));
absorbance_samples = cell(limits2(1)/nreps,limits2(2));
absorbance_std = cell(limits2(1)/nreps,limits2(2));
for i = 1:limits2(1)/nreps
    for j = 1:limits2(2)
        tempNum = limits2(1)/cic;
%         tempVals = zeros(length(absorbance_corrected{1,1}),tempNum);
        tempVals = zeros(length(absorbance_corrected{i,j}),tempNum);
        for k = 1:tempNum
            tempVals(:,k) = absorbance_corrected{tempNum*(i-1)+k,j};
%             disp(k);
        end
        % wrong cases to be deleted
        if delWrong == 1
%             if setup.enzymeName == 'gapdhr' % old implementation
            if contains(setup.enzymeName,"ald") %ph7.81,DF1,rII had an error
                if((i==1)&&(j==7))
                    tempVals = [tempVals(:,1), tempVals(:,2)];
                end
            elseif contains(setup.enzymeName,"gapdhr") %ph7.81,DF1,rII had an error
                if((i==2)&&(j==20))
                    tempVals = [tempVals(:,1), tempVals(:,3)];
                elseif((i==2)&&(j==12))
                    tempVals = [tempVals(:,2), tempVals(:,3)];
                end
            elseif contains(setup.enzymeName,"pgi") %ph6.60,DF1,rI had an error
                if((i==1)&&(j==12))
                    tempVals = [tempVals(:,2), tempVals(:,3)];
                end
            elseif contains(setup.enzymeName,"pfk") %ph7.81,DF16,rIII had an error
                if((i==1)&&(j==22))
                    tempVals = [tempVals(:,1), tempVals(:,2)];
                end
            elseif contains(setup.enzymeName,"pgm") %ph6.19,DF8,rIII had an error
                if((i==1)&&(j==3))
                    tempVals = [tempVals(:,1), tempVals(:,2)];
                end
            elseif contains(setup.enzymeName,"pdc") %multiple deviations, from the initial concentration but similar profile
                if((i==2)&&(j==24))
                    tempVals = [tempVals(:,1), tempVals(:,2)];
                elseif((i==2)&&(j==22))
                    tempVals = [tempVals(:,1), tempVals(:,2)];
                elseif((i==2)&&(j==2))
                    tempVals = [tempVals(:,1), tempVals(:,2)];
                elseif((i==2)&&(j==4))
                    tempVals = [tempVals(:,1), tempVals(:,3)];
                elseif((i==2)&&(j==20))
                    tempVals = [tempVals(:,2), tempVals(:,3)];
                elseif((i==2)&&(j==19))
                    tempVals = [tempVals(:,1), tempVals(:,3)];
                elseif((i==2)&&(j==18))
                    tempVals = [tempVals(:,1), tempVals(:,3)];
                end
            elseif contains(setup.enzymeName,"eno") %multiple deviations, from the initial concentration but similar profile
                if((i==2)&&(j==12))
                    tempVals = [tempVals(:,1), tempVals(:,2)];
                end
            end
        end
        absorbance_mean{i,j} = mean(tempVals,2);
        absorbance_samples{i,j} = tempVals;
        absorbance_std{i,j} = std(tempVals,0,2);
%         disp(i);disp(j);
    end
end

% concentration mean and standard deviation
concentration_mean = cell(limits2(1)/3,limits2(2));
concentration_std = concentration_mean;
for i = 1:limits2(1)/nreps
    for j = 1:limits2(2)
        concentration_mean{i,j} = absorbance_mean{i,j}*extinction_coefficient;
        concentration_std{i,j} = absorbance_std{i,j}*extinction_coefficient;
    end
end

% reaction rates
reaction_rate = cell(limits2(1)/nreps,limits2(2));
for i = 1:limits2(1)/nreps
    for j = 1:limits2(2)
        time_profile = time{excelSheet(i*nreps,j)/2};
        reaction_rate{i,j} = gradient(concentration_mean{i,j})./gradient(time_profile);
    end
end

% reducing size of pH and dilution array
pH_corrected = zeros(limits2(1)/nreps, limits2(2));
Dilution_corrected = pH_corrected;
excelSheet_corrected = pH_corrected;
for i = 1:limits2(1)/nreps
    for j = 1:limits2(2)
        pH_corrected(i,j) = pH{nreps*i,j};
        Dilution_corrected(i,j) = Dilution(nreps*i,j);
        excelSheet_corrected(i,j) = excelSheet(nreps*i,j);
    end
end

% time in cell
time_specific = cell(limits2(1)/3,limits2(2));
for i = 1:limits2(1)/nreps
    for j = 1:limits2(2)
        time_specific{i,j} = time{excelSheet_corrected(i,j)/2};
    end
end

% overal units and protein
unitsTime = 's';
unitsAbs = 'UA';
unitsConc = 'mM';
unitsRates = 'UA s^{-1}';
concProtein = setup.concProtein;
unitsProtein = setup.unitsProtein;
extraDF = setup.extraDF;


% (5) Visualization
tempText = [setup.enzymeName, ' Step 5 of 6: Visualization'];
waitbar(5/6,f,tempText);
pause(1)
pHtested = setup.pHtested;
fullpHarray = setup.fullpHarray;
pHIDs = find(pHtested==1);
pHarray = fullpHarray(pHIDs);
% pHarray = unique(pH_corrected);
enzName = setup.enzymeName;

numpHtested = nnz(pHtested);

if setup.plotOutput == 1
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:numpHtested
        subplot(3,4,i)
        pHval = pHarray(i);
        tempID = find(pH_corrected==pHval);
% % % %         plotConc = zeros(length(concentration_mean{1}),length(tempID));
        plotConc = zeros(length(concentration_mean{tempID(1)}),length(tempID));
        plotStd = plotConc;
        for j = 1:length(tempID)
            plotTime = time{excelSheet_corrected(tempID(j))/2};
            plotConc(:,j) = concentration_mean{tempID(j)};
            plotStd(:,j) = concentration_std{tempID(j)};
% % % %             plotConc(:,j) = absorbance_mean{tempID(j)};
% % % %             plotStd(:,j) = absorbance_std{tempID(j)};
            errorbar(plotTime,plotConc(:,j),plotStd(:,j))
            
            hold on
        end
        % ylim
        if setup.caseStudyENO == 1
            ylim([0 1.5])
        else
            ylim([0 0.15])
        end
        if i == numpHtested
            legend(num2str(Dilution_corrected(tempID)),'location','south','Orientation','horizontal')
        end
        title(erase(sprintf('pH = %d', pHval),"0000e+00"))
    end
    suptitleName = ['Enzyme ', enzName, ': Concentration profile'];
    suptitle(suptitleName);

    cols = ['m','y','r','b'];
    h2 = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:numpHtested
        subplot(3,4,i)
        pHval = pHarray(i);
        tempID = find(pH_corrected==pHval);
% % % % % %         plotConc = zeros(length(concentration_mean{1}),length(tempID));
% %         plotConc = zeros(length(concentration_mean{tempID(1)}),length(tempID));
% %         plotStd = plotConc;

        if setup.caseStudyTPI == 1 % toubleshoot for tpi case
            tempID = tempID(1:4);
        end
        
        for j = 1:length(tempID)
            TData = time{excelSheet_corrected(tempID(j))/2};
            YData = absorbance_samples{tempID(j)};
            stdshade(extinction_coefficient*YData',0.1,cols(j),TData');
%             tempFig = stdshade(extinction_coefficient*YData',0.1,cols(j),TData');
            hold on
        end
        % ylim
        if setup.caseStudyENO == 1
            ylim([0 1.5])
        else
            ylim([0 0.15])
        end
        if i == numpHtested
            % legend development is still not properly finished, right now
            % stdshade cannot be written down so that an output is
            % generated. I cannot get the figures out then.
            
%             tempfig = gcf;
%             set(get(get(tempfig(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%             set(get(get(tempfig(3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%             set(get(get(tempfig(5),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%             set(get(get(tempfig(7),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%             legend('del',num2str(Dilution_corrected(tempID(1))),num2str(Dilution_corrected(tempID(2))),num2str(Dilution_corrected(tempID(3))),num2str(Dilution_corrected(tempID(4))))
            
%             Dilution_corrected(tempID)
% lfact = legend('del','2','del','4','del','6','del','8');
% lfact = legend('del',num2str(Dilution_corrected(tempID(1))),'del',num2str(Dilution_corrected(tempID(2))),'del',num2str(Dilution_corrected(tempID(3))),'del',num2str(Dilution_corrected(tempID(4))));
            if (length(Dilution_corrected(tempID)) == 4)
                legend('del',num2str(Dilution_corrected(tempID(1))),'del',num2str(Dilution_corrected(tempID(2))),'del',num2str(Dilution_corrected(tempID(3))),'del',num2str(Dilution_corrected(tempID(4))));
            elseif (length(Dilution_corrected(tempID)) == 3)
                legend('del',num2str(Dilution_corrected(tempID(1))),'del',num2str(Dilution_corrected(tempID(2))),'del',num2str(Dilution_corrected(tempID(3))));
            elseif (length(Dilution_corrected(tempID)) == 2)
                legend('del',num2str(Dilution_corrected(tempID(1))),'del',num2str(Dilution_corrected(tempID(2))));
            elseif (length(Dilution_corrected(tempID)) == 1)
                legend('del',num2str(Dilution_corrected(tempID(1))));
            end

%               delete(hicot(ismember(get(hicot, 'String'), {'data10','data11','data12'})));
%               delete(lfact(ismember(get(lfact, 'String'), {'del','4'})));


% %             h = zeros(1,3);
% %             legend(h(2:3)); % Only display last two legend titles
% %             legend(h(2:3),'This two','This three');
% %             
% %             legend
%             legend(num2str(Dilution_corrected(tempID)),'location','south','Orientation','horizontal')
        end
        title(erase(sprintf('pH = %d', pHval),"0000e+00"))
    end
    suptitleName = ['Enzyme ', enzName, ': Concentration profile'];
    suptitle(suptitleName);    
    
    h3 = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:numpHtested
        subplot(3,4,i)
        pHval = pHarray(i);
        tempID = find(pH_corrected==pHval);
% % % %         plotRate = zeros(length(concentration_mean{1}),length(tempID));
        plotRate = zeros(length(concentration_mean{tempID(1)}),length(tempID));
        for j = 1:length(tempID)
            plotTime = time{excelSheet_corrected(tempID(j))/2};
            plotRate(:,j) = reaction_rate{tempID(j)};
            plot(plotTime,plotRate(:,j))
            hold on
        end
        if setup.caseStudyPDC == 1
        elseif setup.caseStudyENO == 1
            ylim([0 0.0015])
        else
            ylim([-0.002 0])
        end
        if i == numpHtested
            legend(num2str(Dilution_corrected(tempID)),'location','south','Orientation','horizontal')
        end
        title(erase(sprintf('pH = %d', pHval),"0000e+00"))
    end
    suptitleName = ['Enzyme ', enzName, ': Reaction rate'];
    suptitle(suptitleName);
end

% (6) Saving output variables and figures
tempText = [setup.enzymeName, ' Step 6 of 6: Saving output'];
waitbar(6/6,f,tempText);
pause(1)
%rawData
output.rawData.excelSheet = excelSheet;
output.rawData.coordinates = coordinates;
output.rawData.pH = pH;
output.rawData.dilution = Dilution;
output.rawData.absorbance_raw = absorbance_raw;
output.rawData.absorbance_background = absorbance_background;
output.rawData.time = time;
output.rawData.time_background = time_background;
output.rawData.absorbance_corrected = absorbance_corrected;
%treatedData
output.treatedData.excelSheet_corrected = excelSheet_corrected;
output.treatedData.pH_corrected = pH_corrected;
output.treatedData.dilution_corrected = Dilution_corrected;
output.treatedData.absorbance_mean = absorbance_mean;
output.treatedData.absorbance_samples = absorbance_samples;
output.treatedData.absorbance_std = absorbance_std;
output.treatedData.concentration_mean = concentration_mean;
output.treatedData.concentration_std = concentration_std;
output.treatedData.time = time_specific;
output.treatedData.reaction_rate = reaction_rate;
output.treatedData.unitsTime = unitsTime;
output.treatedData.unitsAbsorbance = unitsAbs;
output.treatedData.unitsConcentration = unitsConc;
output.treatedData.unitsRates = unitsRates;
output.treatedData.concProtein = concProtein;
output.treatedData.unitsProtein = unitsProtein;
output.treatedData.protDF = extraDF;   

if setup.saveOutput == 1
    tempFolderName = ['data/processed_data/',enzName];
    if(~exist(tempFolderName,'dir'))
%         mkdir results/gapdhr
        mkdir('data\processed_data',enzName)
    end
    saveName = ['data\processed_data\', enzName,'\', enzName, '_output.mat'];
    saveFigName1 = ['data\processed_data\', enzName,'\', enzName, '_concentrations.fig'];
    saveFigName2 = ['data\processed_data\', enzName,'\', enzName, '_concentrations_stdshade.fig'];
    saveFigName3 = ['data\processed_data\', enzName,'\', enzName, '_fluxes.fig'];

    %saving
    save(saveName,'output');
    if setup.plotOutput == 1
        savefig(h1,saveFigName1);
        savefig(h2,saveFigName2);
        savefig(h3,saveFigName3);
    end
end

close(f);

end

% %%
% figure
% 
% subplot(1,6,1)
% plot(time{3},absorbance_corrected{4,11},time{3},absorbance_corrected{5,11},time{3},absorbance_corrected{6,11})
% legend('r.I','r.II','r.III')
% 
% subplot(1,6,2)
% plot(time{3},absorbance_corrected{4,12},time{3},absorbance_corrected{5,12},time{3},absorbance_corrected{6,12})
% legend('r.I','r.II','r.III')
% 
% subplot(1,6,3)
% plot(time{6},absorbance_corrected{4,13},time{6},absorbance_corrected{5,13},time{6},absorbance_corrected{6,13})
% legend('r.I','r.II','r.III')
% 
% subplot(1,6,4)
% plot(time{6},absorbance_corrected{4,14},time{6},absorbance_corrected{5,14},time{6},absorbance_corrected{6,14})
% legend('r.I','r.II','r.III')
% 
% subplot(1,6,5)
% plot(time{6},absorbance_corrected{4,15},time{6},absorbance_corrected{5,15},time{6},absorbance_corrected{6,15})
% legend('r.I','r.II','r.III')
% 
% subplot(1,6,6)
% plot(time{6},absorbance_corrected{4,16},time{6},absorbance_corrected{5,16},time{6},absorbance_corrected{6,16})
% legend('r.I','r.II','r.III')

% % %%
% % 
% figure
% 
% subplot(2,2,1)
% plot(time{3},absorbance_raw{4,11},time{3},absorbance_raw{5,11},time{3},absorbance_raw{6,11})
% legend('r.I','r.II','r.III')
% 
% subplot(2,2,2)
% plot(time{3},absorbance_background_interp{4,11},time{3},absorbance_background_interp{5,11},time{3},absorbance_background_interp{6,11})
% legend('r.I','r.II','r.III')
% 
% subplot(2,2,3)
% plot(time{3},background{4,11},time{3},background{5,11},time{3},background{6,11})
% legend('r.I','r.II','r.III')
% 
% subplot(2,2,4)
% plot(time{3},absorbance_corrected{4,11},time{3},absorbance_corrected{5,11},time{3},absorbance_corrected{6,11})
% legend('r.I','r.II','r.III')
% 
% 

% figure
% 
% subplot(2,2,1)
% plot(time_background{(excelSheet2(i,j)+1)/2})
% title('time, background')
% 
% subplot(2,2,2)
% plot(absorbance_background{i,j})
% title('abs, background')
% 
% subplot(2,2,3)
% plot(time{excelSheet(i,j)/2})
% title('time')
% 
% subplot(2,2,4)
% plot(time{4},absorbance_background_interp{i,j})
% title('absorbance, background, interp')
% 
% figure
% plot(output.rawData.time{4},output.rawData.absorbance_corrected{4,12})
% hold on
% plot(output.rawData.time{4},output.rawData.absorbance_corrected{5,12})
% hold on
% plot(output.rawData.time{4},output.rawData.absorbance_corrected{6,12}) % this is the wrong replicate
% legend('r.I','r.Ii','r.III')
%
% figure % better nothing to change in this one
% plot(output.rawData.time{4},output.rawData.absorbance_corrected{1,24})
% hold on
% plot(output.rawData.time{4},output.rawData.absorbance_corrected{2,24})
% hold on
% plot(output.rawData.time{4},output.rawData.absorbance_corrected{3,24}) % this is the wrong replicate
% legend('r.I','r.Ii','r.III')
