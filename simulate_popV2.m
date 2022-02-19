function [Tpop_states,Tpop_propensities,T_rxns]=simulate_popV2(G,Tin_pop_states,t_run) %,display_rxn_on,realtime_print_sim_progress)

display_rxn_on=false; %for now

%% Take this section out when you make this a function
%ngenes=size(G.Nodes,1);
%ncells=size(;
%Tpop_states=array2table(false([ncells,ngenes])) %start at zero state

%%

%assert(islogical(display_rxn_on))

initialize_size_of_rxns_list=100000; %100,000 = initialize list length of reaction events (to make time_list, rxn_ID_list, cell_of_rxn_list)

if istable(Tin_pop_states)
    
%     if ~all(strcmp(G.Nodes.Name',Tin_pop_states.Properties.VariableNames))
%     error('input table Tin_pop_states should have same VariableNames as the G.Nodes table')
%     end

    col_indices_of_Tin_pop_states_are_genes=ismember(Tin_pop_states.Properties.VariableNames,G.Nodes.Name');
    
    if sum(col_indices_of_Tin_pop_states_are_genes==1)~=height(G.Nodes)
       error('did not find the names of G.Nodes.Name in the Tin_pop_states VariableNames') 
    end
    
    % make sure that the column variables are in the same order as G.Nodes.Name
    if ~all(strcmp(Tin_pop_states.Properties.VariableNames(col_indices_of_Tin_pop_states_are_genes),G.Nodes.Name'))
        error('make sure that the column variables of Tin_pop_states in the same order as G.Nodes.Name')
    end
    
    pop_states=Tin_pop_states{:,col_indices_of_Tin_pop_states_are_genes}; %there can be more columns of Tin_pop_states than just genes
    
    %pop_states=Tin_pop_states{:,:}; %convert if table to logical table
else
    pop_states=Tin_pop_states;
end

%% make initial propensities matrix
assert(islogical(pop_states))

ncells=size(pop_states,1);
ngenes=size(G.Nodes,1);

assert(ngenes==size(pop_states,2))

nreactions=ngenes*2; %production and decay for both
pop_propensities=nan(ncells,nreactions);

% call update_state_and_propensities in 'Update_propensities_for_genes'
% mode for all genes
for icell=1:ncells   
        state=pop_states(icell,:);
        propensities=pop_propensities(icell,:);

        % update propensities for all genes
        [state,propensities]=update_state_and_propensities(G.Nodes,state,propensities,'Update_propensities_for_genes',1:ngenes);
        
        pop_states(icell,:)=state;
        pop_propensities(icell,:)=propensities;
end
%Tpop_states{:,:}=pop_states;



%%
t=0; %initialize
rxn_num=0; %initialize
time_list=nan(initialize_size_of_rxns_list,1);
rxn_ID_list=nan(initialize_size_of_rxns_list,1);
cell_of_rxn_list=nan(initialize_size_of_rxns_list,1);

% simulate_pop
while true
    
    if rxn_num==initialize_size_of_rxns_list
        % need to make record-keeping lists even larger - go twice as long
        time_list=[time_list;nan(initialize_size_of_rxns_list,1)];
        rxn_ID_list=[rxn_ID_list;nan(initialize_size_of_rxns_list,1)];
        cell_of_rxn_list=[cell_of_rxn_list;nan(initialize_size_of_rxns_list,1)];
        
        initialize_size_of_rxns_list=2*initialize_size_of_rxns_list;
    end
    
    propensities_sum=sum(sum(pop_propensities));
    
    % waiting time at which this reaction is to occur
    tau=-1/propensities_sum * log(rand());
    t=t+tau;
    
    if t>t_run
        break
    end
    
    % number of this reaction
    rxn_num=rxn_num+1;
    
    time_list(rxn_num)=t;
    
    % reaction that occurs at this time
    rand_reaction_num=rand();%from 0 to 1
    
    [cell_index,rxn_ID]=choose_reaction(pop_propensities,rand_reaction_num);
    
    % record this reaction and the cell it's occuring in
    rxn_ID_list(rxn_num)=rxn_ID;
    cell_of_rxn_list(rxn_num)=cell_index;
    
    %update state and propensities
    old_state=pop_states(cell_index,:);
    old_propensities=pop_propensities(cell_index,:);
    
    [new_state,new_propensities]=update_state_and_propensities(G.Nodes,old_state,old_propensities,'Run_rxn',rxn_ID);
    %%%% [new_state,new_propensities]=updateV1_4binary_states(old_state,old_propensities,rxn_ID);
    
    if display_rxn_on
        fprintf('----- rxn_num= %i, cell_index=%i, rxn_ID = %i, tau= %3.3f , t= %3.3f  ----------------------------\n',rxn_num,cell_index,rxn_ID,tau,t)
        fprintf('',rxn_ID);
        if rxn_ID<=4,Rtype='pro';else Rtype='deg';end
        fprintf('               '); fprintf(sprintf(repmat('       ',[1,rem(rxn_ID-1,4)]))); fprintf('%s\n',Rtype),
        display_state_and_propensities(old_state,old_propensities)
        display_state_and_propensities(new_state,new_propensities)
        fprintf('\n')
    end
    pop_states(cell_index,:)=new_state;
    pop_propensities(cell_index,:)=new_propensities;
    
end

%% truncate lists to right size
rxn_ID_list=rxn_ID_list(~isnan(rxn_ID_list));
cell_of_rxn_list=cell_of_rxn_list(~isnan(cell_of_rxn_list));
time_list=time_list(~isnan(time_list));
%%
T_rxns=table([1:rxn_num]',time_list,cell_of_rxn_list,rxn_ID_list);
T_rxns.Properties.VariableNames={'rxn_num','t','cell_of_rxn','rxn_ID'};

Tpop_states=Tin_pop_states; %initialize for sake of non-gene columns;
Tpop_states(:,col_indices_of_Tin_pop_states_are_genes)=array2table(pop_states);%,'VariableNames',Tin_pop_states.Properties.VariableNames);

assert(all(strcmp(Tpop_states.Properties.VariableNames(col_indices_of_Tin_pop_states_are_genes),G.Nodes.Name')))

Tpop_propensities=array2table(pop_propensities); %first ngenes columns are for production, last ngenes columns is for decay
Tpop_propensities.Properties.VariableNames=[strcat(G.Nodes.Name','_pro'),strcat(G.Nodes.Name','_dec')];



end