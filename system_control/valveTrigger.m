function valveTrigger(~, ~, a, tStart, ax1)

global step_timepoint

fprintf(a,'first_on/');           % Turn on the valve drivers

set(ax1,'Color',[1 0 0 0.2])

step_timepoint = [step_timepoint toc(tStart)];

pause(10)

fprintf(a,'first_off/');           % Turn off the valve drivers

set(ax1,'Color','w')

step_timepoint = [step_timepoint toc(tStart)];

end