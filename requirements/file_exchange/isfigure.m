function OK = isfigure(h)
if strcmp(get(h,'type'),'figure')
  OK = 1;
else
  OK = 0;
end