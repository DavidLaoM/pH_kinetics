fig = figure('Name','Numbers');
t = uitable(fig,'Data',[2,4,6,8;1,3,5,7]); 
saveas(gcf,'Numbers','pdf')