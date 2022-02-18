% Produce pdf from the plots generated
clear,

% % create the pdfs
% export_fig
% print

% % concatenate pdfs
% append_pdfs

% print

openfig('ald_concentrations_stdshade.fig');
print('ald_concentrations_stdshade.pdf','-dpdf','-fillpage')

openfig('gapdh_concentrations_stdshade.fig');
print('gapdh_concentrations_stdshade.pdf','-dpdf','-fillpage')

openfig('gapdhr_concentrations_stdshade.fig');
print('gapdhr_concentrations_stdshade.pdf','-dpdf','-fillpage')

openfig('hxk_concentrations_stdshade.fig');
print('hxk_concentrations_stdshade.pdf','-dpdf','-fillpage')

openfig('pdc_concentrations_stdshade.fig');
print('pdc_concentrations_stdshade.pdf','-dpdf','-fillpage')

openfig('pfk_concentrations_stdshade.fig');
print('pfk_concentrations_stdshade.pdf','-dpdf','-fillpage')

openfig('pgi_concentrations_stdshade.fig');
print('pgi_concentrations_stdshade.pdf','-dpdf','-fillpage')

openfig('pgm_concentrations_stdshade.fig');
print('pgm_concentrations_stdshade.pdf','-dpdf','-fillpage')

openfig('pyk_concentrations_stdshade.fig');
print('pyk_concentrations_stdshade.pdf','-dpdf','-fillpage')

% h=gcf;
% set(h,'Position',[50 50 1200 800]);
% set(h,'PaperOrientation','landscape');
% print('tpi_concentrations.pdf','-dpdf','-fillpage')
openfig('tpi_concentrations.fig');
print('tpi_concentrations.pdf','-dpdf','-fillpage')


%% memoryDump
% % ald
% load('ald_output.mat');
% setup.caseStudyALD = 1;
% setup.caseStudyENO = 0;
% setup.caseStudyGAPDH = 0;
% setup.caseStudyGAPDHr = 0;
% setup.caseStudyHXK = 0;
% setup.caseStudyPDC = 0;
% setup.caseStudyPFK = 0;
% setup.caseStudyPGI = 0;
% setup.caseStudyPGM = 0;
% setup.caseStudyPYK = 0;
% setup.caseStudyTPI = 0;
% 
% %select data
% selectSetup_pH;
% % pHdata
% pHtested = setup.pHtested;
% fullpHarray = setup.fullpHarray;
% pHIDs = find(pHtested==1);
% pHarray = fullpHarray(pHIDs);
% % pHarray = unique(pH_corrected);
% enzName = setup.enzymeName;
% numpHtested = nnz(pHtested);
% 
% % recall data
% 
% %rawData
% excelSheet = output.rawData.excelSheet;
% coordinates = output.rawData.coordinates;
% pH = output.rawData.pH;
%     pH2 = cell2mat(pH);
% Dilution = output.rawData.dilution;
% absorbance_raw = output.rawData.absorbance_raw;
% absorbance_background = output.rawData.absorbance_background;
% time = output.rawData.time;
% time_background = output.rawData.time_background;
% absorbance_corrected = output.rawData.absorbance_corrected;
% %treatedData
% excelSheet_corrected = output.treatedData.excelSheet_corrected;
% pH_corrected = output.treatedData.pH_corrected;
% Dilution_corrected = output.treatedData.dilution_corrected;
% absorbance_mean = output.treatedData.absorbance_mean;
% absorbance_std = output.treatedData.absorbance_std;
% concentration_mean = output.treatedData.concentration_mean;
% concentration_std = output.treatedData.concentration_std;
% time_specific = output.treatedData.time;
% reaction_rate = output.treatedData.reaction_rate;
% unitsTime = output.treatedData.unitsTime;
% unitsAbs = output.treatedData.unitsAbsorbance;
% unitsConc = output.treatedData.unitsConcentration;
% unitsRates = output.treatedData.unitsRates;
% concProtein = output.treatedData.concProtein;
% unitsProtein = output.treatedData.unitsProtein;
% extraDF = output.treatedData.protDF;   
% 
% cols = ['m','y','r','b'];
% 
% h2 = figure('units','normalized','outerposition',[0 0 1 1]);
% for i = 1:numpHtested
%     subplot(3,4,i)
%     pHval = pHarray(i);
% % %     tempID = find(pH_corrected==pHval);
%     tempID = find(pH2==pHval);
% % %     plotConc = zeros(length(concentration_mean{tempID(1)}),length(tempID));
% % %     plotStd = plotConc;
%     for j = 1:length(tempID)
% % %         plotTime = time{excelSheet_corrected(tempID(j))/2};
% % %         plotConc(:,j) = concentration_mean{tempID(j)};
% % %         plotStd(:,j) = concentration_std{tempID(j)};
% % % %         errorbar(plotTime,plotConc(:,j),plotStd(:,j))
%         TData = time(excelSheet(j/2));
%         YData = 1;
%         
% %         plot(Tdata,Ydata,'color','b')
%         stdshade(YData',0.1,cols(j),TData')
%         
% 
%         hold on
%     end
%     ylim([0 0.15])
%     if i == numpHtested
%         legend(num2str(Dilution_corrected(tempID)),'location','south','Orientation','horizontal')
%     end
%     title(erase(sprintf('pH = %d', pHval),"0000e+00"))
% end
% suptitleName = ['Enzyme ', enzName, ': Concentration profile'];
% suptitle(suptitleName);
