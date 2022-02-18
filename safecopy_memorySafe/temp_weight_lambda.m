% weight_lambda


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
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfunction,x_temp,lb,ub,options,data,setup);
    t = toc;
    [raw_error] = costfunction(xres,data,setup);
    setup.selectedLambda = 1;
    [error] = costfunction(xres,data,setup);
    % seting in output arrays
    array_xres{i} = xres;
    array_eData{i} = error(1:end-4);
    array_eParams{i} = error(end-3:end);
    array_resnorm{i} = resnorm;
    array_residual{i} = residual;
    array_Jacobian{i} = Jacobian;
    array_raw_error{i} = raw_error;
end