function [expData] = concatenateImpData(data)
%UNTITLED this file gets all the imported data and writes it down in a 
% sigle data structure.
%   Detailed explanation goes here

fn = fieldnames(data);
fn_lowercase = cell(length(fn),1);
expData_pre = cell(length(fn),1);
fn_cell = struct2cell(data);
% for each term in the structure 'data'
for k=1:numel(fn)
%     enzName = fn{k};
    fn_lowercase{k} = lower(fn{k});
    expData_pre{k} = fn_cell{k};    
%     saveName
%     saveLoc
%     saveV
end
expData = cell2struct(expData_pre, fn_lowercase, 1);
save('data/processed_data/expData.mat','expData');


% % ALD
% load('ald_output.mat','output');
% expData.ald = output;
% % ENO
% load('eno_output.mat','output');
% expData.eno = output;
% % GAPDH
% load('gapdh_output.mat','output');
% expData.gapdh = output;
% % GAPDHr
% load('gapdhr_output.mat','output');
% expData.gapdhr = output;
% % HXK
% load('hxk_output.mat','output');
% expData.hxk = output;
% % PDC
% load('pdc_output.mat','output');
% expData.pdc = output;
% % PFK
% load('pfk_output.mat','output');
% expData.pfk = output;
% % PGI
% load('pgi_output.mat','output');
% expData.pgi = output;
% % PGM
% load('pgm_output.mat','output');
% expData.pgm = output;
% % PYK
% load('pyk_output.mat','output');
% expData.pyk = output;
% % TPI
% load('tpi_output.mat','output');
% expData.tpi = output;
% 
% save('results/expData.mat','expData');

end

