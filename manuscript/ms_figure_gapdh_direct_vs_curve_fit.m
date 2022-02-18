% c_royalBlue = [65	105	225]/255; % royalblue
% c_midnightblue = [25	25	112]/255; % midnightblue

%% 1 Recalling only GAPFH fwd figure

enzymeName = {'gapdhr'}; %6 gapdh_fwd

for i = 1:length(enzymeName)
    % close to call by the same numbers
    close all
    % load the figure + rescale
    loadName_train = [enzymeName{i}, '_trainData_fit_metabolites_reg.fig'];
    fh1 = openfig(loadName_train);
    set(gcf, 'units','normalized','outerposition',[0 0 0.5 1]);
    fh1_children = get(fh1,'Children');
    for j = 1:length(fh1_children)
        if fh1_children(j).Tag == "suptitle"
            delete(fh1_children(j));
        else
            fh1_children(j).Children(1).Position(2) = fh1_children(j).YLim(2)*0.85;
            fh1_children(j).Children(1).Position(1) = fh1_children(j).XLim(2)*0.25;
        end
    end
    fh1.Position = [fh1.Position(1) fh1.Position(2) 0.9*fh1.Position(3) 0.9*fh1.Position(4) ];
end


%% 2 selecting only the pH7.06 and plotting in a separate cell
% %% identification of case 
for i = [1,3:10] % it's number 6
    % 
    disp(i);
    disp(fh1_children(i).Position);
end
% %% copy and relocate figure
f2 = figure();
ax2 = copyobj(fh1_children(6),f2);
f2.Position = [100 100 800 800];
ax2.Position = [0.25 0.40 0.50 0.50];
% %% selecting and deleting extra
ax2_children = get(ax2,'Children');
delete(ax2_children(9))
delete(ax2_children(8))
delete(ax2_children(7))
delete(ax2_children(6))
delete(ax2_children(3))
delete(ax2_children(2))
delete(ax2_children(1))


% %% 3 adding also the direct method plot
hold on
temp_x = [ax2_children(4).XData(1), 300];
temp_slope = (ax2_children(4).YData(2) - ax2_children(4).YData(1)) / ax2_children(4).XData(2);
temp_y = [ax2_children(4).YData(1), ax2_children(4).YData(1) + temp_slope * 300];
% % %% plot
% plot(temp_x, temp_y, 'r-', 'LineWidth',2)
% %% plot and add color
c_royalBlue = [65	105	225]/255; % royalblue
c_midnightblue = [25	25	112]/255; % midnightblue
plot(temp_x, temp_y, '-', 'LineWidth',2,'color',c_royalBlue)
ax2_children(5).Color = c_midnightblue;

% %% 4 adding text labels and color
% specific
% text(300/50, 0.15/25, 'V_{GAPDHR} = V_{GAPDHR.max}')
% % long_equation = [{'v_{GAPDHfwd = (Vmf .* (GAP .* NAD - BPG .* NADH ./ Keq)./(Kgap .* Knad))./...
% %     ((1 + NAD ./ Knad + NADH ./ Knadh) .* (1 + BPG ./...
% %     Kbpg + GAP ./ Kgap));
% 
%     
% long_equation = 'v_{GAPDHR} = ((V_{GAPDHR.max} .* (BPG .* NADH ./ K_{eq} - GAP .* NAD))';
% % './(Kbpg .* Knadh))./...
% %     ((1 + NAD ./ Knad + NADH ./ Knadh) .* (1 + BPG ./...
% %     Kbpg + GAP ./ Kgap));
% 
% long_equation = '^{(BPG .* NADH ./ K_{eq} - GAP .* NAD)}/_{b}';
% 
% text(300/10*5, 0.15/10*5, long_equation)


% Build a string that contains the Latex expression
eqtext = '$$V_{GAPDHR}=V_{GAPDHR}';
eqtext = [eqtext '\left ( {    BPG .* NADH ./ Keq - GAP .* NAD  \over   (Kbpg .* Knadh) .*((1 + NAD ./ Knad + NADH ./ Knadh) .* (1 + BPG ./ Kbpg + GAP ./ Kgap))}   \right)^n$$'];
eqtext2 = '$$V_{GAPDHR}=V_{GAPDHR}$$';

% Add the string containing the Latex expression to the plot
text(125, -0.15/3, eqtext2, 'Interpreter', 'Latex', 'FontSize', 10, 'Color', 'k')
text(-125, -0.15/2, eqtext, 'Interpreter', 'Latex', 'FontSize', 10, 'Color', 'k')

% general
xlabel('Time (s)')
ylabel('NADH concentration (mM)')
ax2.FontSize = 12;

% color white
set(f2,'color','w')

% legend
% delete(hL1)
% 
hL1 = legend(f2.Children, 'Curve fitting simulations', ...
                            'Experimental data2',...
                            'Direct fitting simulations');
% hL1.Orientation = 'horizontal';
hL1.Box = 'off';
% hL1.FontSize = 11;
% hL1.Position = [0.03    0.04    0.3960    0.0340];

%% saving
if setup.saveOutput == 1
    % 
    savefig(f2,'1appendix_GAPDRH_dir_curve')
    % 
    set(f2,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters 1appendix_GAPDRH_dir_curve
    print -dpng -painters 1appendix_GAPDRH_dir_curve
end
