enzymeName = setup.enzymeName;
% parameter estimation
array_xres = cell(1,length(lambdalist));
array_eData = cell(1,length(lambdalist));
array_eParams = cell(1,length(lambdalist));
    array_resnorm = cell(1,length(lambdalist));
    array_residual = cell(1,length(lambdalist));
    array_Jacobian = cell(1,length(lambdalist));
    array_raw_error = cell(1,length(lambdalist));
for i = 1:length(lambdalist)
    fprintf('pEst for lambda=%d\n',lambdalist(i));
    setup.selectedLambda = lambdalist(i);
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t = toc;
    [raw_error] = optfun(xres,data,setup);
    setup.selectedLambda = 1;
    [error] = optfun(xres,data,setup);
    % seting in output arrays
    array_xres{i} = xres;
    switch enzymeName
        case 'hxk'
% %             array_eData{i} = error(1:end-16);
% %             array_eParams{i} = error(end-15:end);
%             array_eData{i} = error(1:end-16);
%             array_eParams{i} = error(end-15:end-12);
            array_eData{i} = error(1:end-4);
            array_eParams{i} = error(end-3:end);
        case 'pgi'
% %             array_eData{i} = error(1:end-18);
% %             array_eParams{i} = error(end-17:end);
%             array_eData{i} = error(1:end-14);
%             array_eParams{i} = error(end-13:end-12);
            array_eData{i} = error(1:end-2);
            array_eParams{i} = error(end-1:end);
        case 'pfk'
% %             array_eData{i} = error(1:end-25);
% %             array_eParams{i} = error(end-24:end);
%             array_eData{i} = error(1:end-25);
%             array_eParams{i} = error(end-24:end-12);
            array_eData{i} = error(1:end-13);
            array_eParams{i} = error(end-12:end);
        case 'ald'
% %             array_eData{i} = error(1:end-13);
% %             array_eParams{i} = error(end-12:end);
%             array_eData{i} = error(1:end-13);
%             array_eParams{i} = error(end-12:end-10);
            array_eData{i} = error(1:end-3);
            array_eParams{i} = error(end-2:end);
        case 'tpi'
% %             array_eData{i} = error(1:end-9);
% %             array_eParams{i} = error(end-8:end);
%             array_eData{i} = error(1:end-9);
%             array_eParams{i} = error(end-8:end-7);
            array_eData{i} = error(1:end-2);
            array_eParams{i} = error(end-1:end);
        case 'gapdh'
% %             array_eData{i} = error(1:end-16);
% %             array_eParams{i} = error(end-15:end);
%             array_eData{i} = error(1:end-16);
%             array_eParams{i} = error(end-15:end-12);
%             array_eData{i} = error(1:end-4);
%             array_eParams{i} = error(end-3:end);
            switch caseKm
                case 'pH_independent'
                    array_eData{i} = error(1:end-4);
                    array_eParams{i} = error(end-3:end);
                case 'pH_dependent'
                    array_eData{i} = error(1:end-48);
                    array_eParams{i} = error(end-47:end);
            end
            
        case 'gapdhr'
% %             array_eData{i} = error(1:end-14);
% %             array_eParams{i} = error(end-13:end);
%             array_eData{i} = error(1:end-14);
%             array_eParams{i} = error(end-13:end-10);
%             array_eData{i} = error(1:end-4);
%             array_eParams{i} = error(end-3:end);
            switch caseKm
                case 'pH_independent'
                    array_eData{i} = error(1:end-4);
                    array_eParams{i} = error(end-3:end);
                case 'pH_dependent'
                    array_eData{i} = error(1:end-40);
                    array_eParams{i} = error(end-39:end);
            end
        
        case 'eno'
% %             array_eData{i} = error(1:end-14);
% %             array_eParams{i} = error(end-13:end);
%             array_eData{i} = error(1:end-14);
%             array_eParams{i} = error(end-13:end-12);
            array_eData{i} = error(1:end-2);
            array_eParams{i} = error(end-1:end);
        case 'pgm'
% %             array_eData{i} = error(1:end-14);
% %             array_eParams{i} = error(end-13:end);
%             array_eData{i} = error(1:end-14);
%             array_eParams{i} = error(end-13:end-12);
            array_eData{i} = error(1:end-2);
            array_eParams{i} = error(end-1:end);
        case 'pyk'
% %             array_eData{i} = error(1:end-18);
% %             array_eParams{i} = error(end-17:end);
%             array_eData{i} = error(1:end-18);
%             array_eParams{i} = error(end-17:end-12);
            array_eData{i} = error(1:end-6);
            array_eParams{i} = error(end-5:end);
        case 'pdc'
% %             array_eData{i} = error(1:end-14);
% %             array_eParams{i} = error(end-13:end);
%             array_eData{i} = error(1:end-14);
%             array_eParams{i} = error(end-13:end-12);
            array_eData{i} = error(1:end-2);
            array_eParams{i} = error(end-1:end);
        case 'gapdh_FWD_REV'
            array_eData{i} = error(1:end-49);
            array_eParams{i} = error(end-47:end);
        otherwise
            disp('Warning: eData and eParams has not been selected');
    end
        array_resnorm{i} = resnorm;
        array_residual{i} = residual;
        array_Jacobian{i} = Jacobian;
        array_raw_error{i} = raw_error;
end