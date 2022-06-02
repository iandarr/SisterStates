% Gillespie to simulate gene networks and infer direct biomolecular
% interactions based on sister-cell differential analysis
% 
% gene exprssion states are binary (on=1, off=0)

%% simulation inputs
clear

simulation_file='sim1_parallel_posFeedback.mat';
rng(1)
P.ncells=10000;
P.t_run_steady_state=20; %run for t_run time units;
display_rxn_on=false;
realtime_print_sim_progress=20; %seconds (actual, real time seconds) report back on simluation progress and expected simulation time after this many

% %% Create gene network G (mixture of edge relationships)
% s = {'U',   'X',    'Y',    'U',    'P1',   'P2', 'R'};
% t = {'X',   'Y',    'Z',    'P1',   'P2',   'P3', 'X'};
% G = digraph(s,t); %create digraph
% %                       U-X      X-Y    Y-Z     U-P1    P1-P2   P2-P3   R-X      Ualt-X
% G.Edges.rate_affecting={'pro', 'pro',  'pro',  'pro',  'dec',  'pro',  'dec'}';
% G.Edges.funct_type=    {'mult', 'mult','add',  'add',  'mult', 'add',  'add'}';
% G.Edges.funct_mag=     [10,    10,      1,      1,      10,     -1,      -2]';

% %% Create gene network G (linear)
% s = {'g1',   'g2',    'g3',    'g4',    'g5',   'g6', 'g7', 'g8', 'g9'};
% t = {'g2',   'g3',    'g4',    'g5',   'g6',   'g7', 'g8', 'g9', 'g10'};
% G = digraph(s,t); %create digraph
% %                       U-X      X-Y    Y-Z     U-P1    P1-P2   P2-P3   R-X      Ualt-X
% G.Edges.rate_affecting=repmat({'dec'},[numedges(G),1]);%{'pro', 'pro',  'pro',  'pro',  'pro',  'pro',  'pro'}';
% G.Edges.funct_type=    repmat({'mult'},[numedges(G),1]);% {'mult', 'mult','mult',  'mult',  'mult', 'mult',  'mult'}';
% G.Edges.funct_mag=     repmat(0.2,[numedges(G),1]);%[10,    10,      10,      10,      10,     10,      10]';

% %% Create gene network G (parallel)
% s = {'U',   'x1',    'x2',    'x3',    'x4',    'x5',   'x6',    'U',    'p1',    'p2',    'p3',    'p4',   'p5',   'p6'};
% t = {'x1',   'x2',    'x3',    'x4',   'x5',    'x6',   'x7',    'p1',    'p2',    'p3',    'p4',   'p5',   'p6',    'p7'};
% G = digraph(s,t); %create digraph
% %                       U-X      X-Y    Y-Z     U-P1    P1-P2   P2-P3   R-X      Ualt-X
% G.Edges.rate_affecting=repmat({'dec'},[numedges(G),1]);%{'pro', 'pro',  'pro',  'pro',  'pro',  'pro',  'pro'}';
% G.Edges.funct_type=    repmat({'mult'},[numedges(G),1]);% {'mult', 'mult','mult',  'mult',  'mult', 'mult',  'mult'}';
% G.Edges.funct_mag=     repmat(0.2,[numedges(G),1]);%[10,    10,      10,      10,      10,     10,      10]';

%% Create gene network G (parallel, with postive feedback)
s = {'U',   'x1',    'x2',    'x3',    'x4',    'x5',   'x6',    'U',    'p1',    'p2',    'p3',    'p4',   'p5',   'p6'};
t = {'x1',   'x2',    'x3',    'x4',   'x5',    'x6',   'x7',    'p1',    'p2',    'p3',    'p4',   'p5',   'p6',    'p7'};
s(end+1)={'x4'};
t(end+1)={'x3'};
G = digraph(s,t); %create digraph
%                       U-X      X-Y    Y-Z     U-P1    P1-P2   P2-P3   R-X      Ualt-X
G.Edges.rate_affecting=repmat({'dec'},[numedges(G),1]);%{'pro', 'pro',  'pro',  'pro',  'pro',  'pro',  'pro'}';
G.Edges.funct_type=    repmat({'mult'},[numedges(G),1]);% {'mult', 'mult','mult',  'mult',  'mult', 'mult',  'mult'}';
G.Edges.funct_mag=     repmat(0.2,[numedges(G),1]);%[10,    10,      10,      10,      10,     10,      10]';

% x4-->x3 affects production (current code can't handle 2 edges affecting a given rate (pro or dec)
idx=all([strcmp(G.Edges.EndNodes(:,1),'x4'),strcmp(G.Edges.EndNodes(:,2),'x3')],2);
G.Edges.rate_affecting(idx)={'pro'};
G.Edges.funct_type(idx)={'mult'};
G.Edges.funct_mag(idx)=4;
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
G.Nodes.NatProRate(1)=0.2; %first node
G.Nodes.NatDecRate(1)=0.2; %first node

G.Edges.relationship_label=join([G.Edges.rate_affecting,G.Edges.funct_type,cellstr(num2str(G.Edges.funct_mag))]); %just for the graph, not used elsewhere
G.Edges.relationship_label_short=replace(G.Edges.relationship_label,{'pro','dec','mult','add',' '},{'P','D','*','+',''});
% now populate the G.Nodes table with columns that indicate predecessors and
% matrices (lookup tables) of how the propensities are updated based the predecessors 'preds'
% G.Nodes.proPreds, G.Nodes.proMatrix, G.Nodes.decPreds, G.Nodes.decMatrix
G=make_node_rate_matrices(G);

% plot the network graph
clf
plot(G,'EdgeLabel',G.Edges.relationship_label_short)

%% Initialize population by simulating it to a 'steady state'
ngenes=size(G.Nodes,1);
Tinit_pop_states=array2table(false([P.ncells,ngenes])); %start at all zero states
Tinit_pop_states.Properties.VariableNames=G.Nodes.Name';
% SIMULATE POPULATION
fprintf("simulating network to a 'steady state'...")
[Tss_pop_states,Tss_pop_propensities,T_rxns]=simulate_popV2(G,Tinit_pop_states,P.t_run_steady_state);%,display_rxn_on,realtime_print_sim_progress);
fprintf("done\n")
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
P.ncells=height(Tdiv_pop_states);
% %% add some labels to the population
Tdiv_cellLabels=table([1:P.ncells]',floor([1:.5:P.ncells/2+.5])','VariableNames',{'cellInd','parentInd'});
Tdiv_pop_states=[Tdiv_cellLabels,Tdiv_pop_states]; %add cell labels

%% Now run the population in short incremental periods which are shorter than the timescale of a single reaction
% since production reactions occur at rate 1 in our updateV1_4binary
% function, want to run it less than this.
P.timepoint_interval=0.05;
P.maxSortTime=1;
P.time_to_run=2;
tnow=0;

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
%alpha=0.05; % alpha for confidence interval on plot

P.sortGenes=G.Nodes.Name; % use all genes to sort cells by
%GroupingGenes=G.Nodes.Name; % use all genes to group cells by

S=struct();
iS=0;
P.randSeedStart=0;
%Meas(iMeas).data=DifferentialAnalysisBinofit(Tdiv_pop_states,GroupingGenes);
%Meas(iMeas).time=tnow;
while tnow<=P.time_to_run
    fprintf('tnow=%f of %f\n',tnow,P.time_to_run)
    
    % make a hypothetical 'sort' (don't really need to do this for all
    % times though
    if tnow<=P.maxSortTime
    iS=iS+1;
    for iSortGene=1:length(P.sortGenes)
        sortGene=P.sortGenes{iSortGene};
        S(iS).(sortGene).idxHi=Tdiv_pop_states.(sortGene);
        S(iS).(sortGene).dtSort=tnow;
    end
    end
    
    %% make measurements on sorted cells.
    %iMeas=iMeas+1;
    % first decide on past sort times to take measurements for
    iStoMeasureList=1:length(S); % all of them
    % now take measurements (Both sister-naive and sister-control measurements)
    for  iStoMeasure=iStoMeasureList
        sortGeneList=fields(S(iStoMeasure))';
        for iSortGene=1:length(sortGeneList)
            sortGene=sortGeneList{iSortGene};
            idxHiAtSortTime=S(iStoMeasure).(sortGene).idxHi;
            P.randSeedStart=P.randSeedStart+1;
        	Meas=DifferentialAnalysisBinofit2(Tdiv_pop_states,idxHiAtSortTime,P.randSeedStart);
            if ~isfield(S(iStoMeasure).(sortGene),'meas')
                iMeas=1;
            else
                iMeas=length(S(iStoMeasure).(sortGene).meas)+1;
            end
            S(iStoMeasure).(sortGene).meas(iMeas).data=Meas;
            % record times
            S(iStoMeasure).(sortGene).meas(iMeas).time=tnow;
            dtSort=S(iStoMeasure).(sortGene).dtSort;
            S(iStoMeasure).(sortGene).meas(iMeas).dtMeas=tnow-dtSort;
        end
    end
    
    % simulate the population forward
    if tnow<=P.time_to_run
        [Tdiv_pop_states,Tdiv_pop_propensities,Tdiv_rxns]=simulate_popV2(G,Tdiv_pop_states,P.timepoint_interval);
        tnow=tnow+P.timepoint_interval;% update the time
    end
end
%% save/load the simulation (sort S, parameters P, initial population states Tinit_pop_states, and gene network G)
save(simulation_file,'P','G','Tinit_pop_states','S')

%% load and view a simulation
simulation_file='sim1_parallel_posFeedback.mat';
load(simulation_file)
ViewSisterStates(S,G)





