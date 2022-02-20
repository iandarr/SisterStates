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

% %% Create gene network G
% s = {'U',   'X',    'Y',    'U',    'P1',   'P2', 'R'};
% t = {'X',   'Y',    'Z',    'P1',   'P2',   'P3', 'X'};
% G = digraph(s,t); %create digraph
% %                       U-X      X-Y    Y-Z     U-P1    P1-P2   P2-P3   R-X      Ualt-X
% G.Edges.rate_affecting={'pro', 'pro',  'pro',  'pro',  'dec',  'pro',  'dec'}';
% G.Edges.funct_type=    {'mult', 'mult','add',  'add',  'mult', 'add',  'add'}';
% G.Edges.funct_mag=     [10,    10,      1,      1,      10,     -1,      -2]';

%% Create gene network G (linear)
s = {'g1',   'g2',    'g3',    'g4',    'g5',   'g6', 'g7', 'g8', 'g9'};
t = {'g2',   'g3',    'g4',    'g5',   'g6',   'g7', 'g8', 'g9', 'g10'};
G = digraph(s,t); %create digraph
%                       U-X      X-Y    Y-Z     U-P1    P1-P2   P2-P3   R-X      Ualt-X
G.Edges.rate_affecting=repmat({'pro'},[numedges(G),1]);%{'pro', 'pro',  'pro',  'pro',  'pro',  'pro',  'pro'}';
G.Edges.funct_type=    repmat({'mult'},[numedges(G),1]);% {'mult', 'mult','mult',  'mult',  'mult', 'mult',  'mult'}';
G.Edges.funct_mag=     repmat(10,[numedges(G),1]);%[10,    10,      10,      10,      10,     10,      10]';

% %% Create gene network G (parallel)
% s = {'U',   'x1',    'x2',    'x3',    'x4',    'x5',   'x6',    'U',    'p1',    'p2',    'p3',    'p4',   'p5',   'p6'};
% t = {'x1',   'x2',    'x3',    'x4',   'x5',    'x6',   'x7',    'p1',    'p2',    'p3',    'p4',   'p5',   'p6',    'p7'};
% G = digraph(s,t); %create digraph
% %                       U-X      X-Y    Y-Z     U-P1    P1-P2   P2-P3   R-X      Ualt-X
% G.Edges.rate_affecting=repmat({'pro'},[numedges(G),1]);%{'pro', 'pro',  'pro',  'pro',  'pro',  'pro',  'pro'}';
% G.Edges.funct_type=    repmat({'mult'},[numedges(G),1]);% {'mult', 'mult','mult',  'mult',  'mult', 'mult',  'mult'}';
% G.Edges.funct_mag=     repmat(10,[numedges(G),1]);%[10,    10,      10,      10,      10,     10,      10]';

%%
nnodes=numnodes(G);
nedges=numedges(G);
G.Edges.funct_priority=nan(nedges,1);
G.Edges.funct_priority(1)=1;
G.Nodes.Properties.RowNames=G.Nodes.Name; %creates a referencable row name

% indicate the basal production and decay rates
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

%% plot the network graph
clf
plot(G,'EdgeLabel',G.Edges.relationship_label_short)

%% Initialize population by simulating it to a 'steady state'
ngenes=size(G.Nodes,1);
Tinit_pop_states=array2table(false([ncells,ngenes])); %start at all zero states
Tinit_pop_states.Properties.VariableNames=G.Nodes.Name';
%% SIMULATE POPULATION2
[Tss_pop_states,Tss_pop_propensities,T_rxns]=simulate_popV2(G,Tinit_pop_states,t_run_steady_state);%,display_rxn_on,realtime_print_sim_progress);
% or load steady state from file
%clear
%load mySS_50

%% Add table labels and divide population
% 1 division with perfect heritability of state

%Tdiv_pop_states=array2table(nan(2*height(Tss_pop_states),width(Tss_pop_states)));
Tdiv_pop_states=table();
%Tdiv_pop_states.Properties.VariableNames=Tss_pop_states.Properties.VariableNames;
Tdiv_pop_states(1:2:2*height(Tss_pop_states),:)=Tss_pop_states; % odd rows
Tdiv_pop_states(2:2:2*height(Tss_pop_states),:)=Tss_pop_states; % even rows
ncells=height(Tdiv_pop_states);
% %% add some labels to the population
Tdiv_cellLabels=table([1:ncells]',floor([1:.5:ncells/2+.5])','VariableNames',{'cellInd','parentInd'});
Tdiv_pop_states=[Tdiv_cellLabels,Tdiv_pop_states]; %add cell labels

%% Now run the population in short incremental periods which are shorter than the timescale of a single reaction
% since production reactions occur at rate 1 in our updateV1_4binary
% function, want to run it less than this.
timepoint_interval=.5;
time_to_run=10;
tnow=0;
rng(1)
% %% Now simulate forward in time. At longer periods of time, we will
% simulate the results from taking single-cell measurements and calculating
% one of the following:
% (a) sister-naive calculation: for cells that are high for a given gene (Eg. cells where g1=1), subtract
%       their other genes' average states (Eg. g2, g3, g4,...) from the g1=0
%       cells' states.
% (b) sister-control calculation: for cells that are high for a given gene
%       (Eg. cells where g1=1) where their sister cells are low got that gene
%       (g1=0), subtract their other genes' average states from their sister's
%       states, and find the average of those differences.

% timepoints and monitoring for how states change
MeasureSisterNaive=struct();
MeasureSisterControl=struct();
alpha=0.05; % alpha for confidence interval on plot

iMeas=0;
while tnow<=time_to_run
    [Tdiv_pop_states,Tdiv_pop_propensities,Tdiv_rxns]=simulate_popV2(G,Tdiv_pop_states,timepoint_interval);
    tnow=tnow+timepoint_interval;% update the time

    % take hypothetical measurements
    iMeas=iMeas+1;
    
    % take sister naive and sister control measurements
    GroupingGenes=G.Nodes.Name; % use all genes to group cells by
    Meas(iMeas).data=DifferentialAnalysisBinofit(Tdiv_pop_states,GroupingGenes);
    Meas(iMeas).time=tnow;
end

%%
tutorialApp1(Meas)




