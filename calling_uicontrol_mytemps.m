%%
subH=subplot(2,2,[1 3]);
get(gca,'position')

%%
figure
uipanel_handle = uipanel(gcf);
%%
mytemps(pos)

%%
figure
plot([1 2],[3 4])
ax=gca;
aspect = get(ax,'PlotBoxAspectRatio');
% Change axes Units property (this only works with non-normalized units)
set(ax,'Units','pixels');
% Get and update the axes width to match the plot aspect ratio.
pos = get(ax,'Position');
pos(3) = aspect(1)/aspect(2)*pos(4);
set(ax,'Position',pos);

%%
f = gcf;
c = uicontrol(f,'Style','popupmenu');

