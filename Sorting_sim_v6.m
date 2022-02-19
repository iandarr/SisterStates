% Gillespie to simulate gene networks in a population of cells
% gene exprssion states are binary (on=1, off=0)
% v6 WANT to update (still untouched)
%   sorting structure S is S(1), S(2)... and has .settings, .data, and .summary
% v5 updated:
%   to plot mean and error bars (based on confidence interval from binofit function)
%   and put the plotting into a function, population_plot.m
% v4 updated:
%   include a directed graph gene network input 'G' (copied from graph_test2.m)
%   but a given gene can only have up to one incoming edge that
%   affects its production 'pro' rate and one incoming edge that affects its
%   decay 'dec' rate.

% Wish list of updates:
%       -faster simulations. Consider simulate_pop.m to accept
%       T_propensities (can fill in first time if it's empty)
%       -put DT into S, and make S cell array of structures S(1), S(2) each
%        and with Flag1, Flag2, Flag3 (can sort after Flag2 and Flag3)
%       with sort settings information. like this: S(1).data.g2=[table]
%       -ploting can accept multiple sorts (with interactive menu) and multiple
%       genes to sort based on (with interactive menu)
%       input/output files
%       - digraph is plotted horizontally, and red arrows for repressive
%       - print periodic simulation progress
%       - natural prod/decay rates are shown in table in figure
%% simulation inputs
clear
rng(1)
ncells=5;
t_run_steady_state=50; %run for t_run time units;
display_rxn_on=false;
realtime_print_sim_progress=20; %seconds (actual, real time seconds) report back on simluation progress and expected simulation time after this many

%% Create gene network G
s = {'U',   'X',    'Y',    'U',    'P1',   'P2', 'R'};
t = {'X',   'Y',    'Z',    'P1',   'P2',   'P3', 'X'};
G = digraph(s,t); %create digraph
%                       U-X      X-Y    Y-Z     U-P1    P1-P2   P2-P3   R-X      Ualt-X
G.Edges.rate_affecting={'pro', 'pro',  'pro',  'pro',  'dec',  'pro',  'dec'}';
G.Edges.funct_type=    {'mult', 'mult','add',  'add',  'mult', 'add',  'add'}';
G.Edges.funct_mag=     [10,    10,      1,      1,      10,     -1,      -2]';

%% Create gene network G (linear)
s = {'g1',   'g2',    'g3',    'g4',    'g5',   'g6', 'g7', 'g8', 'g9'};
t = {'g2',   'g3',    'g4',    'g5',   'g6',   'g7', 'g8', 'g9', 'g10'};
G = digraph(s,t); %create digraph
%                       U-X      X-Y    Y-Z     U-P1    P1-P2   P2-P3   R-X      Ualt-X
G.Edges.rate_affecting=repmat({'pro'},[numedges(G),1]);%{'pro', 'pro',  'pro',  'pro',  'pro',  'pro',  'pro'}';
G.Edges.funct_type=    repmat({'mult'},[numedges(G),1]);% {'mult', 'mult','mult',  'mult',  'mult', 'mult',  'mult'}';
G.Edges.funct_mag=     repmat(10,[numedges(G),1]);%[10,    10,      10,      10,      10,     10,      10]';

%% Create gene network G (parallel)
s = {'U',   'x1',    'x2',    'x3',    'x4',    'x5',   'x6',    'U',    'p1',    'p2',    'p3',    'p4',   'p5',   'p6'};
t = {'x1',   'x2',    'x3',    'x4',   'x5',    'x6',   'x7',    'p1',    'p2',    'p3',    'p4',   'p5',   'p6',    'p7'};
G = digraph(s,t); %create digraph
%                       U-X      X-Y    Y-Z     U-P1    P1-P2   P2-P3   R-X      Ualt-X
G.Edges.rate_affecting=repmat({'pro'},[numedges(G),1]);%{'pro', 'pro',  'pro',  'pro',  'pro',  'pro',  'pro'}';
G.Edges.funct_type=    repmat({'mult'},[numedges(G),1]);% {'mult', 'mult','mult',  'mult',  'mult', 'mult',  'mult'}';
G.Edges.funct_mag=     repmat(10,[numedges(G),1]);%[10,    10,      10,      10,      10,     10,      10]';

%%
nnodes=numnodes(G);
nedges=numedges(G);
G.Edges.funct_priority=nan(nedges,1);
G.Edges.funct_priority(1)=1;
G.Nodes.Properties.RowNames=G.Nodes.Name; %creates a referencable row name

% indicate the production and decay rates
norm_ProRate=1;
norm_DecRate=1;
G.Nodes.NatProRate=ones(nnodes,1)*norm_ProRate; %apply general rates
G.Nodes.NatDecRate=ones(nnodes,1)*norm_DecRate; %apply general rates
G.Nodes.NatProRate(1)=1; %first node
G.Nodes.NatDecRate(1)=1; %first node

G.Edges.relationship_label=join([G.Edges.rate_affecting,G.Edges.funct_type,cellstr(num2str(G.Edges.funct_mag))]); %just for the graph, not used elsewhere
G.Edges.relationship_label_short=replace(G.Edges.relationship_label,{'pro','dec','mult','add',' '},{'P','D','*','+',''});
% now populate the G.Nodes table with columns that indicate predecessors and
% matrices (lookup tables) of how the propensities are updated based the predecessors 'preds'
% G.Nodes.proPreds, G.Nodes.proMatrix, G.Nodes.decPreds, G.Nodes.decMatrix

G=make_node_rate_matrices(G);

%%
clf
plot(G,'EdgeLabel',G.Edges.relationship_label_short)

%% Initialize population at a 'steady state'

ngenes=size(G.Nodes,1);
Tinit_pop_states=array2table(false([ncells,ngenes])); %start at zero state
Tinit_pop_states.Properties.VariableNames=G.Nodes.Name';

%% SIMULATE POPULATION2

[Tss_pop_states,Tss_pop_propensities,T_rxns]=simulate_popV2(G,Tinit_pop_states,t_run_steady_state);%,display_rxn_on,realtime_print_sim_progress);

%% or load from file
%clear
%load mySS_50

%% Add table labels and divide population
% %% simulate 1 division with perfect heritability of state

%Tdiv_pop_states=array2table(nan(2*height(Tss_pop_states),width(Tss_pop_states)));
Tdiv_pop_states=table();
%Tdiv_pop_states.Properties.VariableNames=Tss_pop_states.Properties.VariableNames;
Tdiv_pop_states(1:2:2*height(Tss_pop_states),:)=Tss_pop_states;
Tdiv_pop_states(2:2:2*height(Tss_pop_states),:)=Tss_pop_states;
ncells=height(Tdiv_pop_states);
% %% add some labels to the population
Tdiv_cellLabels=table([1:ncells]',floor([1:.5:ncells/2+.5])','VariableNames',{'cellInd','parentInd'});
Tdiv_pop_states=[Tdiv_cellLabels,Tdiv_pop_states]; %add cell labels

%%
sortGeneNames=G.Nodes.Name;% sort based on all of them

S=struct(); 
S(1).settings.sortGeneNames=G.Nodes.Name;
S(1).settings.Flag1='lo-lo';
S(1).settings.Flag2='hi-lo';
S(1).settings.Flag3=[];
S(1).settings.SortTrigger=[]; %always Flag2 as written
S(1).settings.WaitTimeAfterSortTrigger=1;
S(1).settings.AbortEvent='AbortSortIfChangedWhileWaiting';
% %S(1).data=[];
% S(2).settings.sortGeneNames=G.Nodes.Name;
% S(2).settings.Flag1='lo-lo';
% S(2).settings.Flag2='hi-lo';
% S(2).settings.Flag3=[];
% S(2).settings.SortTrigger=[]; %always Flag2 as written
% S(2).settings.WaitTimeAfterSortTrigger=0;
% S(2).settings.AbortEvent='AbortSortIfChangedWhileWaiting';

%%
tnow=100;
S=sister_sort_v2(S,Tdiv_pop_states,tnow)

%% Now run the population in short incremental periods which are shorter than the timescale of a single reaction
% since production reactions occur at rate 1 in our updateV1_4binary
% function, want to run it less than this.
timepoint_interval=.5;
time_to_run=10;
tnow=0;
rng(1)
% %% Now simulate forward in time, 'scanning the plate' at different
% timepoints and monitoring for how states change
while tnow<=time_to_run
    [Tdiv_pop_states,Tdiv_pop_propensities,Tdiv_rxns]=simulate_popV2(G,Tdiv_pop_states,timepoint_interval);
    tnow=tnow+timepoint_interval;% update the time
    %% sister sorting
    S=sister_sort_v2(S,Tdiv_pop_states,tnow);%,sortType,twait,sortGeneNames);
    
end

%% Get a Distribution Table (DT)
% create structure Dsum that summarizes the average expression states for
% poppulations. S=sorted (sisters); B=Bulk, at end of simulation:
% Sall, Shi, Slo, Ball, Bhi, Blo
alpha=0.05; % alpha for confidence interval on plot

DT=struct();
nSortGenes=length(sortGeneNames);
for sortGeneInd=1:nSortGenes
% create empty elements in DT
sortGeneName=sortGeneNames{sortGeneInd};
DT.(sortGeneName)=array2table(nan(6+12,1+nSortGenes)); % 1 column for nCells; 12 columns for confidence interval
DT.(sortGeneName).Properties.VariableNames=[{'nCells'};G.Nodes.Name];
DT.(sortGeneName).Properties.RowNames={'Sall','Shi', 'Slo', 'Ball', 'Bhi', 'Blo'...
    'Sall_botCI','Sall_topCI','Shi_botCI','Shi_topCI','Slo_botCI','Slo_topCI','Ball_botCI','Ball_topCI','Bhi_botCI','Bhi_topCI','Blo_botCI','Blo_topCI'}; % added 17mar19
% add in bottom and top of confidence intervals of proportion
%DT.(sortGeneName).Properties.RowNames=[DT.(sortGeneName).Properties.RowNames,{'Shi_botCI','Shi_topCI','Slo_botCI','Slo_topCI','Bhi_botCI','Bhi_topCI','Blo_botCI','Blo_topCI'}];

%% first look at sister sorted structure S
isSorted=~isnan(S.(sortGeneName).SortTime);
isSortedAndFlaggedHi=all([strcmp(S.(sortGeneName).FlagType,'hi'),~isnan(S.(sortGeneName).SortTime)],2);
isSortedAndFlaggedLo=all([strcmp(S.(sortGeneName).FlagType,'lo'),~isnan(S.(sortGeneName).SortTime)],2);

DT.(sortGeneName){{'Sall','Shi','Slo'},'nCells'}=[sum(isSorted),sum(isSortedAndFlaggedHi),sum(isSortedAndFlaggedLo)]';

%testMat=S.(sortGeneName){isSorted,G.Nodes.Name};
%DT.(sortGeneName){'Sall_OLD',G.Nodes.Name}=mean(S.(sortGeneName){isSorted,G.Nodes.Name});
%DT.(sortGeneName){'Shi_OLD',G.Nodes.Name}=mean(S.(sortGeneName){isSortedAndFlaggedHi,G.Nodes.Name});
%DT.(sortGeneName){'Slo_OLD',G.Nodes.Name}=mean(S.(sortGeneName){isSortedAndFlaggedLo,G.Nodes.Name});

% Sorted all cells
[proportions,CIs]=binofit(sum(S.(sortGeneName){isSorted,G.Nodes.Name})',sum(isSorted),alpha);
DT.(sortGeneName){'Sall',G.Nodes.Name}=proportions';
DT.(sortGeneName){'Sall_botCI',G.Nodes.Name}=CIs(:,1)';
DT.(sortGeneName){'Sall_topCI',G.Nodes.Name}=CIs(:,2)';

% Sorted hi cells (with 95% CI)
[proportions,CIs]=binofit(sum(S.(sortGeneName){isSortedAndFlaggedHi,G.Nodes.Name})',sum(isSortedAndFlaggedHi),alpha);
DT.(sortGeneName){'Shi',G.Nodes.Name}=proportions';
DT.(sortGeneName){'Shi_botCI',G.Nodes.Name}=CIs(:,1)';
DT.(sortGeneName){'Shi_topCI',G.Nodes.Name}=CIs(:,2)';

% Sorted lo cells (with 95% CI)
[proportions,CIs]=binofit(sum(S.(sortGeneName){isSortedAndFlaggedLo,G.Nodes.Name})',sum(isSortedAndFlaggedLo),alpha);
DT.(sortGeneName){'Slo',G.Nodes.Name}=proportions';
DT.(sortGeneName){'Slo_botCI',G.Nodes.Name}=CIs(:,1)';
DT.(sortGeneName){'Slo_topCI',G.Nodes.Name}=CIs(:,2)';

%% now look at bulk population (Tdiv_pop_states)

isBulkHi=Tdiv_pop_states.(sortGeneName)==1;
isBulkLo=Tdiv_pop_states.(sortGeneName)==0;

%nCells
DT.(sortGeneName){{'Ball','Bhi','Blo'},'nCells'}=[height(Tdiv_pop_states),sum(isBulkHi),sum(isBulkLo)]';

%DT.(sortGeneName){'Ball_OLD',G.Nodes.Name}=mean(Tdiv_pop_states{:,G.Nodes.Name});
%DT.(sortGeneName){'Bhi',G.Nodes.Name}=mean(Tdiv_pop_states{isBulkHi,G.Nodes.Name});
%DT.(sortGeneName){'Blo',G.Nodes.Name}=mean(Tdiv_pop_states{isBulkLo,G.Nodes.Name});

% Bulk full population
[proportions,CIs]=binofit(sum(Tdiv_pop_states{:,G.Nodes.Name},1)',height(Tdiv_pop_states),alpha);
DT.(sortGeneName){'Ball',G.Nodes.Name}=proportions';
DT.(sortGeneName){'Ball_botCI',G.Nodes.Name}=CIs(:,1)';
DT.(sortGeneName){'Ball_topCI',G.Nodes.Name}=CIs(:,2)';

% Bulk hi cells (with 95% CI)
[proportions,CIs]=binofit(sum(Tdiv_pop_states{isBulkHi,G.Nodes.Name},1)',sum(isBulkHi,1),alpha);
DT.(sortGeneName){'Bhi',G.Nodes.Name}=proportions';
DT.(sortGeneName){'Bhi_botCI',G.Nodes.Name}=CIs(:,1)';
DT.(sortGeneName){'Bhi_topCI',G.Nodes.Name}=CIs(:,2)';

% Bulk lo cells (with 95% CI)
[proportions,CIs]=binofit(sum(Tdiv_pop_states{isBulkLo,G.Nodes.Name},1)',sum(isBulkLo,1),alpha);
DT.(sortGeneName){'Blo',G.Nodes.Name}=proportions';
DT.(sortGeneName){'Blo_botCI',G.Nodes.Name}=CIs(:,1)';
DT.(sortGeneName){'Blo_topCI',G.Nodes.Name}=CIs(:,2)';


end

%% Plot based on this structure, Distribution Tables (DT)

sortGeneToPlot='x7';
clf

%%% plot the Gene Interaction Graph
subplot(2,2,[1 3])
plot(G,'EdgeLabel',G.Edges.relationship_label_short)

%%% plot the sister sorted population
subplot(2,2,2)
%pop_this_plot={'Sall','Shi','Slo'};
%mycolors={'black','green','red'};
pop_this_plot={'Shi','Slo'};
mycolors={'green','red'};
population_plot(G.Nodes.Name,DT,sortGeneToPlot,pop_this_plot,mycolors)
title(sprintf('Sister sort on %s, timepoint interval=%.2g, twait=%.2g, sortType=%s',sortGeneToPlot,timepoint_interval,twait,sortType),'Interpreter', 'none')

%%% plot the bulk sorted population
subplot(2,2,4)
%pop_this_plot={'Sall','Shi','Slo'};
%mycolors={'black','green','red'};
pop_this_plot={'Bhi','Blo'};
mycolors={'green','red'};
population_plot(G.Nodes.Name,DT,sortGeneToPlot,pop_this_plot,mycolors)
title(sprintf('Bulk sort at end of simulation on %s',sortGeneToPlot))
