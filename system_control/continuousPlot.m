function continuousPlot(src, event, ax1, ax2, ax3, ax4)

global times datas

times = [times; event.TimeStamps];
datas = [datas; event.Data];

xlims = [max([0 ceil(event.TimeStamps(end))-10]) ceil(event.TimeStamps(end))];

plot(ax1, times(max([0 end-10*400])+1:end,1), datas(max([0 end-10*400])+1:end,1), 'b')
ax1.XLim = xlims;
ax1.YLim = [0 2];

speedDeg = diff(datas(max([0 end-10*400])+1:end,2));
plot(ax2, times(max([0 end-10*400])+1:end,1), [speedDeg; speedDeg(end)], 'b');
ax2.XLim = xlims;
ax2.YLim = [-2 2];

plot(ax3, times(max([0 end-10*400])+1:end,1), datas(max([0 end-10*400])+1:end,3), 'b')
ax3.XLim = xlims;
% ax3.YLim = [1.5 4];

plot(ax4, times(max([0 end-10*400])+1:end,1), datas(max([0 end-10*400])+1:end,4), 'b')
ax4.XLim = xlims;
% ax4.YLim = [1 3];

end