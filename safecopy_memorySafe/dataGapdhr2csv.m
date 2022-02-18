% %% (0) Setup and export data 2 csv

clear
set_paths_pHstudy;
dbstop if error

% select specific case and recall data
setup.caseStudyALD = 0;
setup.caseStudyENO = 0;
setup.caseStudyGAPDH = 0;
setup.caseStudyGAPDHr = 1;
setup.caseStudyHXK = 0;
setup.caseStudyPDC = 0;
setup.caseStudyPFK = 0;
setup.caseStudyPGI = 0;
setup.caseStudyPGM = 0;
setup.caseStudyPYK = 0;
setup.caseStudyTPI = 0;
selectSetup_pH;
% added
setup.saveOutput = 0;

load('expData.mat','expData');
import_gapdhR = expData.gapdhr;

DFs = setup.DFactorsTotal;
pHtested = setup.pHtested;
numpHtested = nnz(pHtested);
pHs = numpHtested;
blank = zeros(pHs,DFs);
blankCell = cell(pHs,DFs);

% data reorganization
pHTemp = blank';
DFTemp = blank';
abs_meanTemp = blankCell';
abs_stdTemp = blankCell';
conc_meanTemp = blankCell';
conc_stdTemp = blankCell';
timeTemp = blankCell';
RRsTemp = blankCell';

pHarray = unique(import_gapdhR.treatedData.pH_corrected);
for i = 1:numpHtested
    pHval = pHarray(i);
    tempID = find(import_gapdhR.treatedData.pH_corrected==pHval);
    pHTemp(:,i) = import_gapdhR.treatedData.pH_corrected(tempID);
    DFTemp(:,i) = import_gapdhR.treatedData.dilution_corrected(tempID);
    for j = 1:4
        abs_meanTemp{j,i} = import_gapdhR.treatedData.absorbance_mean{tempID(j)};
        abs_stdTemp{j,i} = import_gapdhR.treatedData.absorbance_std{tempID(j)};
        conc_meanTemp{j,i} = import_gapdhR.treatedData.concentration_mean{tempID(j)};
        conc_stdTemp{j,i} = import_gapdhR.treatedData.concentration_std{tempID(j)};
        timeTemp{j,i} = import_gapdhR.treatedData.time{tempID(j)};
        RRsTemp{j,i} = import_gapdhR.treatedData.reaction_rate{tempID(j)};
    end
end

pH = pHTemp';
DF = DFTemp';
abs_mean = abs_meanTemp';
abs_std = abs_stdTemp';

conc_mean = conc_meanTemp';
conc_std = conc_stdTemp';
time = timeTemp';
RRs = RRsTemp';
Vmax = blank';

NADH = blankCell;
Vmax = blank;
for i = 1:(DFs*numpHtested)
    tempDiff = conc_mean{i} - conc_mean{i}(1); % all stoichiometries are 1-to-1.
    NADH{i} = conc_mean{i};
%     RRs2 = RRs';
    Vmax(i) = max(abs(RRs{i}));
end
pH = pHTemp';
DF = DFTemp';
abs_mean = abs_meanTemp';
abs_std = abs_stdTemp';
conc_mean = conc_meanTemp';
conc_std = conc_stdTemp';
time = timeTemp';
RRs = RRsTemp';
clear pHTemp DFTemp abs_meanTemp abs_stdTemp conc_meanTemp conc_stdTemp timeTemp RRsTemp

% save in data
data.pH = pH;
data.DF = DF;
data.abs_mean = abs_mean;
data.abs_std = abs_std;
data.conc_mean = conc_mean;
data.conc_std = conc_std;
data.time = time;
data.RRs = RRs;
data.Vmax = Vmax;
data.chosenVmax = max(max(Vmax));
data.chosenNADini = 0.15;
temp1 = import_gapdhR.rawData.absorbance_corrected{4,4};
temp2 = import_gapdhR.rawData.absorbance_corrected{5,4};
temp3 = import_gapdhR.rawData.absorbance_corrected{6,4};
data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
data.raw.time = import_gapdhR.rawData.time{1};

pHvals = unique(import_gapdhR.treatedData.pH_corrected);
% visualize: check calculations made
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:numpHtested
    subplot(3,4,i)
    for j = 1:DFs
        plot(time{i,j},NADH{i,j},'.-')
        hold on
    end    
    title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
end
suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
suptitle(suptitleName);


%% Exporting to CSV
% writetable
% 'Sheet', 'MySheetName'

T = table(['M';'F';'M'],[45;41;36],...
    {'New York, NY';'San Diego, CA';'Boston, MA'},[true;false;false])
writetable(T,'myData.csv','Delimiter',',','QuoteStrings',true)
type 'myData.csv'

T = table({'M';'M';'F';'F';'F'},[38;43;38;40;49],...
          [71;69;64;67;64],[176;163;131;133;119])
T.Properties.VariableNames = {'Gender' 'Age' 'Height' 'Weight'}


%%
filename = 'gapdhr_data.xlsx';
for i = 1:length(pHarray)
    
    % get data
    DF_1 = 1*ones(size(data.time{i,4}));
    DF_2 = 2*ones(size(data.time{i,3}));
    DF_4 = 4*ones(size(data.time{i,2}));
    DF_8 = 8*ones(size(data.time{i,1}));
    time_1 = data.time{i,4};
    time_2 = data.time{i,3};
    time_4 = data.time{i,2};
    time_8 = data.time{i,1};
    NADHmean_1 = data.conc_mean{i,4};
    NADHmean_2 = data.conc_mean{i,3};
    NADHmean_4 = data.conc_mean{i,2};
    NADHmean_8 = data.conc_mean{i,1};
    NADHstd_1 = data.conc_std{i,4};
    NADHstd_2 = data.conc_std{i,3};
    NADHstd_4 = data.conc_std{i,2};
    NADHstd_8 = data.conc_std{i,1};
    vGAPDH_1 = data.RRs{i,4};
    vGAPDH_2 = data.RRs{i,3};
    vGAPDH_4 = data.RRs{i,2};
    vGAPDH_8 = data.RRs{i,1}; 
    
    % build matrix
    T = table(...
        DF_1, time_1, NADHmean_1, NADHstd_1, vGAPDH_1,...
        DF_2, time_2, NADHmean_2, NADHstd_2, vGAPDH_2,...
        DF_4, time_4, NADHmean_4, NADHstd_4, vGAPDH_4,...
        DF_8, time_8, NADHmean_8, NADHstd_8, vGAPDH_8);
    T.Properties.VariableNames = {...
        'DF1' 'time1' 'NADH_mean1' 'NADH_std1' 'vGAPDH1' ...
        'DF2' 'time2' 'NADH_mean2' 'NADH_std2' 'vGAPDH2' ...
        'DF4' 'time4' 'NADH_mean4' 'NADH_std4' 'vGAPDH4' ...
        'DF8' 'time8' 'NADH_mean8' 'NADH_std8' 'vGAPDH8' ...
        };
    
    % save matrix
    sheetname = erase(sprintf('pH_%d',pHarray(i)),"0000e+00");
    writetable(T,filename,'Sheet',sheetname);
end

