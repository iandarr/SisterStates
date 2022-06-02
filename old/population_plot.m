function linesHandle=population_plot(G_Nodes_Name,DT,sortGeneToPlot,pop_this_plot,mycolors)

assert(length(pop_this_plot)==length(mycolors))

ngenes=length(G_Nodes_Name);

pop_this_plot_botCI=join([pop_this_plot',repmat({'_botCI'},[length(pop_this_plot) 1])],'')';
pop_this_plot_topCI=join([pop_this_plot',repmat({'_topCI'},[length(pop_this_plot) 1])],'')';

Ymeans=DT.(sortGeneToPlot){pop_this_plot,G_Nodes_Name}';
Xdata=ones(ngenes,length(pop_this_plot)).*[1:ngenes]';
Yneg_err=Ymeans - DT.(sortGeneToPlot){pop_this_plot_botCI,G_Nodes_Name}';
Ypos_err=DT.(sortGeneToPlot){pop_this_plot_topCI,G_Nodes_Name}' - Ymeans;
lineHand=errorbar(Xdata,Ymeans,Yneg_err,Ypos_err);

legend_labels=cell(length(pop_this_plot),1);
for iLine=1:numel(lineHand)
    set(lineHand(iLine),'Color',mycolors{iLine})
    legend_labels{iLine}=sprintf('%s (N=%i)',pop_this_plot{iLine},DT.(sortGeneToPlot){pop_this_plot{iLine},'nCells'});
end

legend(legend_labels)
ahandle=gca;
ahandle.XTick=1:length(G_Nodes_Name);
ahandle.XTickLabel=G_Nodes_Name;





end