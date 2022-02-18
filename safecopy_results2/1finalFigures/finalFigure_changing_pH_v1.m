% pH values tested
pH_vals     = [6.19 6.26 6.41 6.60 6.81 7.06 7.29 7.51 7.68 7.81];
% keq manually taken from the equilibrator. All concentrations in [mM]
% units.
kEq.HXK     = [2.74E2 3.03E2 3.79E2 5.16E2 7.45E2 1.2E3 1.9E3 3E3 4.3E3 5.7E3];  %dir+
kEq.PGI     = [3.6E-1 3.6E-1 3.6E-1 3.6E-1 3.61E-1 3.61E-1 3.62E-1 3.63E-1 3.64E-1 3.65E-1];  %dir+
kEq.PFK     = [2.56E1 2.96E1 4.03E1 5.95E1 9.15E1 1.54E2 2.5E2 4.05E2 5.91E2 7.95E2];  %dir+, E.C.num.: EC 2.7.1.56
kEq.FBA     = [1.2E-3 1.1E-3 9.7E-4 8.6E-4 7.9E-4 7.3E-4 7E-4 6.8E-4 6.7E-4 6.6E-4];  %dir+
kEq.TPI     = [1/(8.07) 1/(8.21) 1/(8.46) 1/(8.74) 1/(8.97) 1/(9.16) 1/(9.26) 1/(9.33) 1/(9.36) 1/(9.38)];  %dir-
kEq.GAPDH   = [1.7E-3, 2.3E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1,1.78E-1, 2.45E-1];  %dir+
kEq.PGK     = [1/(7.4E-4), 1/(7.2E-4), 1/(6.9E-4), 1/(6.5E-4), 1/(6.1E-4), 1/(5.7E-4), 1/(5.5E-4), 1/(5.3E-4), 1/(5.2E-4), 1/(5.2E-4)];  %dir-
kEq.ENO     = [5.22 5.22 5.21 5.21 5.2 5.2 5.2 5.19 5.19 5.19];  %dir ?
kEq.PGM     = [1/6.62 1/6.44 1/6.12 1/5.83 1/5.62 1/5.46 1/5.38 1/5.33 1/5.3 1/5.29]; %dir-
kEq.PYK     = [1/(3.1E-6) 1/(3.5E-6) 1/(4.5E-6) 1/(6.5E-6) 1/(9.7E-6) 1/(1.6E-5) 1/(2.5E-5) 1/(4.1E-5) 1/(5.9E-5) 1/(7.9E-5)];  %dir-
kEq.PDC     = [1E4 8.9E3 6.3E3 4.1E3 2.5E3 1.4E3 8.31E2 5.01E2 3.39E2 2.51E2];  %dir+

% Generating the reference values
kEq_ref68 = kEq;
for i = 1:length(pH_vals)
    kEq_ref68.HXK(i) = kEq_ref68.HXK(i)./kEq.HXK(5);     %= [2.74E2 3.03E2 3.79E2 5.16E2 7.45E2 1.2E3 1.9E3 3E3 4.3E3 5.7E3];  %dir+
    kEq_ref68.PGI(i) = kEq_ref68.PGI(i)./kEq.PGI(5);   %= [3.6E-1 3.6E-1 3.6E-1 3.6E-1 3.61E-1 3.61E-1 3.62E-1 3.63E-1 3.64E-1 3.65E-1];  %dir+
    kEq_ref68.PFK(i) = kEq_ref68.PFK(i)./kEq.PFK(5);     %= [2.56E1 2.96E1 4.03E1 5.95E1 9.15E1 1.54E2 2.5E2 4.05E2 5.91E2 7.95E2];  %dir+, E.C.num.: EC 2.7.1.56
    kEq_ref68.FBA(i) = kEq_ref68.FBA(i)./kEq.FBA(5);     %= [1.2E-3 1.1E-3 9.7E-4 8.6E-4 7.9E-4 7.3E-4 7E-4 6.8E-4 6.7E-4 6.6E-4];  %dir+
    kEq_ref68.TPI(i) = kEq_ref68.TPI(i)./kEq.TPI(5);     %= [1/(8.07) 1/(8.21) 1/(8.46) 1/(8.74) 1/(8.97) 1/(9.16) 1/(9.26) 1/(9.33) 1/(9.36) 1/(9.38)];  %dir-
    kEq_ref68.GAPDH(i) = kEq_ref68.GAPDH(i)./kEq.GAPDH(5);   %= [1.7E-3, 2.3E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1,1.78E-1, 2.45E-1];  %dir+
    kEq_ref68.PGK(i) = kEq_ref68.PGK(i)./kEq.PGK(5);     %= [1/(7.4E-4), 1/(7.2E-4), 1/(6.9E-4), 1/(6.5E-4), 1/(6.1E-4), 1/(5.7E-4), 1/(5.5E-4), 1/(5.3E-4), 1/(5.2E-4), 1/(5.2E-4)];  %dir-
    kEq_ref68.ENO(i) = kEq_ref68.ENO(i)./kEq.ENO(5);     %= [5.22 5.22 5.21 5.21 5.2 5.2 5.2 5.19 5.19 5.19];  %dir ?
    kEq_ref68.PGM(i) = kEq_ref68.PGM(i)./kEq.PGM(5);     %= [1/6.62 1/6.44 1/6.12 1/5.83 1/5.62 1/5.46 1/5.38 1/5.33 1/5.3 1/5.29]; %dir-
    kEq_ref68.PYK(i) = kEq_ref68.PYK(i)./kEq.PYK(5);     %= [1/(3.1E-6) 1/(3.5E-6) 1/(4.5E-6) 1/(6.5E-6) 1/(9.7E-6) 1/(1.6E-5) 1/(2.5E-5) 1/(4.1E-5) 1/(5.9E-5) 1/(7.9E-5)];  %dir-
    kEq_ref68.PDC(i) = kEq_ref68.PDC(i)./kEq.PDC(5);     %= [1E4 8.9E3 6.3E3 4.1E3 2.5E3 1.4E3 8.31E2 5.01E2 3.39E2 2.51E2];  %dir+
end

%% Creating the figure
% preexisting
dataPlot = cell(1,length(pH_vals));
dataPlotLeft = dataPlot;
enzymeNames = {'hxk','pgi','pfk','fba','tpi','gapdh','pgk','eno','pgm','pyk','pdc'};
colArr = cool(11);
endVals1 = zeros(11,1);
A1 = struct2cell(kEq_ref68);
A2 = struct2cell(kEq);
for i = 1:11
    endVals1(i) = A1{i}(end);
end
[endVals2,idxs] = sort(endVals1,'ascend');

% figure
f101 = figure(101);

% subplot (right)
subplot(1,3,[2 3])
dataPlot{1} = semilogy(pH_vals,kEq_ref68.HXK,'.-'); hold on, %HXK
dataPlot{2} = semilogy(pH_vals,kEq_ref68.PGI,'.-'); hold on, %PGI
dataPlot{3} = semilogy(pH_vals,kEq_ref68.PFK,'.-'); hold on, %PFK
dataPlot{4} = semilogy(pH_vals,kEq_ref68.FBA,'.-'); hold on, %FBA
dataPlot{5} = semilogy(pH_vals,kEq_ref68.TPI,'.-'); hold on, %TPI
dataPlot{6} = semilogy(pH_vals,kEq_ref68.GAPDH,'.-'); hold on, %GAPDH
dataPlot{7} = semilogy(pH_vals,kEq_ref68.PGK,'.-'); hold on, %PGK
dataPlot{8} = semilogy(pH_vals,kEq_ref68.ENO,'.-'); hold on, %ENO
dataPlot{9} = semilogy(pH_vals,kEq_ref68.PGM,'.-'); hold on, %PGM
dataPlot{10} = semilogy(pH_vals,kEq_ref68.PYK,'.-'); hold on, %PYK
dataPlot{11} = semilogy(pH_vals,kEq_ref68.PDC,'.-'); hold on, %PDC
ylim([1E-2 1E2])
ylabel('Equilibrium constant []')
xlabel('pH []')
for i = 1:length(dataPlot)
    % xpos fixed and loop to get ypos
    xpos = 7.85;
    if i == 1, ypos = dataPlot{i}.YData(end) - 0.25;
    elseif i == 3, ypos = dataPlot{i}.YData(end) + 1.5;
    elseif i == 6, ypos = dataPlot{i}.YData(end) + 4;
    elseif i == 7, ypos = dataPlot{i}.YData(end) + 0.75;
    elseif i == 9, ypos = dataPlot{i}.YData(end) + 0.4;
    elseif i == 2, ypos = dataPlot{i}.YData(end) + 0.15;
    elseif i == 8, ypos = dataPlot{i}.YData(end) - 0.1;
    elseif i == 5, ypos = dataPlot{i}.YData(end) - 0.225;
    elseif i == 4, ypos = dataPlot{i}.YData(end) - 0.25;
    elseif i == 11, ypos = dataPlot{i}.YData(end) - 0.01;
    else ypos = dataPlot{i}.YData(end); end
    text(xpos,ypos,enzymeNames{i});
    dataPlot{idxs(i)}.Color = colArr(i,:);
end

% subplot (left)
subplot(1,3,1)
for i = 1:11
    % xpos fixed and loop to get ypos
    xpos = 1.1;
    ypos = A2{i}(5);
%     if i == 1, ypos = A2{i}(5);
%     elseif i == 3, ypos = A2{i}(5);
%     elseif i == 6, ypos = A2{i}(5);
%     elseif i == 7, ypos = A2{i}(5);
%     elseif i == 9, ypos = A2{i}(5);
%     elseif i == 2, ypos = A2{i}(5);
%     elseif i == 8, ypos = A2{i}(5);
%     elseif i == 5, ypos = A2{i}(5);
%     elseif i == 4, ypos = A2{i}(5);
%     elseif i == 11, ypos = A2{i}(5);
%     else ypos = A2{i}(5); end
    
    dataPlotLeft{i} = semilogy(1,ypos,'o'); hold on
    dataPlotLeft{i}.MarkerFaceColor = colArr(i,:);
    dataPlotLeft{i}.Color = 'black';
    text(xpos,ypos,enzymeNames{i});
%     x2 = [0 2];
%     y2 = [1E-1 1E6];
%     patch(x2,y2,'cool')
%     set(gca,'xtick',[])
end

%%
% t = linspace(0, 2*pi, 100);
% x = cos(t);
% y = sin(t);
% c = x; % colored following x value and current colormap
% 
% figure
% patch(x,y,c)
% % hold on
% % scatter(x,y)
% % hold off
% %%
% x = [2 5; 2 5; 5 8];
% y = [4 0; 8 2; 4 0];
% c = [0 3; 6 4; 4 6];
% % x = [-1 1];
% % y = [-1 1];
% % c = [0 1];
% figure
% patch(x,y,c)
% % colorbar

%%