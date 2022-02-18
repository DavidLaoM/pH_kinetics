function isAxes = IsAxes(ax)
try
    isAxes = strcmp(get(ax, 'type'), 'axes');
catch
    isAxes = false;
end
end