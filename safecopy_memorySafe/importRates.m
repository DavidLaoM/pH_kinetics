% %  importRates.m
% Importing the rection rates that Laura sent on 2020-07-27.
clear,

assayName = {'ald', 'eno', 'gapdh', 'gapdhr', 'hxk', 'pdc', 'pfk', 'pgi', 'pgm', 'pyk', 'tpi'};
lenAssay = length(assayName);
ratesCell = cell(1,lenAssay);

% loop
for i = 1:lenAssay
    meanName = assayName{i}; % selected case
    tempMean = xlsread(meanName); % read mean
    tempStd = xlsread(['std_', meanName]); % read std
    ratesCell{i} = [tempMean(2:end,2:end), tempStd(2:end,2:end)];% locate values in cell
end
% locate values in struct
ratesStruct.ald = ratesCell{1};
ratesStruct.eno = ratesCell{2};
ratesStruct.gapdh = ratesCell{3};
ratesStruct.gapdhr = ratesCell{4};
ratesStruct.hxk = ratesCell{5};
ratesStruct.pdc = ratesCell{6};
ratesStruct.pfk = ratesCell{7};
ratesStruct.pgi = ratesCell{8};
ratesStruct.pgm = ratesCell{9};
ratesStruct.pyk = ratesCell{10};
ratesStruct.tpi = ratesCell{11};

% save results 
save('experimentalRates.mat','assayName','ratesCell','ratesStruct'); % assayName, ratesCell, ratesStruct

%%
% plot to check that everything is fine.
figure,
for i = 1:length(assayName)
    subplot(3,4,i)
    errorbar(ratesCell{i}(:,1),ratesCell{i}(:,2),ratesCell{i}(:,4))
    title(assayName{i})
end



