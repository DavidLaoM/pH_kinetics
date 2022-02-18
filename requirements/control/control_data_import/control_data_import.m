fn = fieldnames(data);
% for each term in the structure 'data'
for k=1:numel(fn)
    % load control
    load(['control_', fn{k} ,'_output.mat']);
    eval(['mat_new = data.',fn{k},';'])
    eval(['mat_control = ',lower(fn{k}),'_output_control;'])
    % check raw data (last in corrected absorbance)
    end_abs_control = mat_control.rawData.absorbance_corrected{end,end}(end);
    end_abs_new = mat_new.rawData.absorbance_corrected{end,end}(end);
    if end_abs_control == end_abs_new
        disp(['Does ',fn{k},' raw data match control safecopy? YES.'])
    else
        disp(['Does ',fn{k},' raw data match control safecopy? NO.'])
    end
    % check imported data (last in pipeline, concentrations)
    end_conc_control = mat_control.treatedData.concentration_mean{end,end}(end);
    end_conc_new = mat_new.treatedData.concentration_mean{end,end}(end);
    if end_conc_control == end_conc_new
        disp(['Does ',fn{k},' processed data match control safecopy? YES.'])
    else
        disp(['Does ',fn{k},' processed data match control safecopy? NO.'])
    end
    
    
end


% % memoryDump
% clear
% 
% enzymeNames = {'ald',...
%                 'eno',...
%                 'gapdh',...
%                 'gapdhr',...
%                 'hxk',...
%                 'pdc',...
%                 'pfk',...
%                 'pgi',...
%                 'pgm',...
%                 'pyk',...
%                 'tpi'};
% for i = 1:length(enzymeNames)
%     loadName = ['control_',enzymeNames{i},'_output.mat'];
%     load(loadName);
%     eval([enzymeNames{i},'_output_control = output; clear output'])
%     save(['control_',enzymeNames{i},'_output.mat'],[enzymeNames{i},'_output_control'])
% end
% clear
% 
% 
% %%
% length(fn)
% 
% % eval([enzymeNames{i},'_output_control = output; clear output'])
% 
% (length(fieldnames(mat_control.rawData))-1) * rand
% 
% 
% figure, histogram()
% 
% (length(fieldnames(mat_control.rawData))-1) * rand+1
% 
% figure, plot(ceil((length(fieldnames(mat_control.rawData))-1) * rand(10000,1)+1))
% 
% arr = round((length(fieldnames(mat_control.rawData))-1) * rand(10000,1)+1);
% figure, plot(ones(10000,1),arr,'o')
% 
% 