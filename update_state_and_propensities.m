function [state_out,propensities_out]=update_state_and_propensities(G_Nodes,state,propensities,varargin)
% [state_out,propensities_out]=update_state_and_propensities(G_Nodes,state,propensities,'Run_rxn',rxn_ID);
% [state_out,propensities_out]=update_state_and_propensities(G_Nodes,state,propensities,'Update_propensities_for_genes',gene_IDs);


% changes to make - most important are succInds, proPreds, decPreds, proMatrix, decMatrix
assert(istable(G_Nodes))

% commenting these assert statements saves 60% of time:
%assert(all(ismember({'Name','NatProRate','NatDecRate','predNames','succNames','proMatrix','decMatrix'},G_Nodes.Properties.VariableNames)))
%assert(all(strcmp(G_Nodes.Properties.RowNames,G_Nodes.Name))) %Names should be RowNames



assert(size(state,1)==1)
assert(size(propensities,1)==1)

ngenes=size(state,2);
nreactions=2*ngenes; %first ngenes columns are production, last ngenes columns are degredation

assert(size(propensities,2)==nreactions)

% % when I do this as cell arrays
% assert(iscell(predInds))
% assert(iscell(succInds))
% assert(size(predInds,1)==ngenes)
% assert(size(succInds,1)==ngenes)
% assert(size(predInds,2)==1)
% assert(size(succInds,2)==1)

if size(varargin,2)==2
   
    update_mode=varargin{1};
    

    switch update_mode
        case 'Run_rxn' %one the one reaction by updating the state and propensities
            rxn_ID=varargin{2};
            if length(rxn_ID)>1, error('can only update one reaction at a time'),end
            
            assert(rem(rxn_ID,1)==0) %whole numbers
            assert(all(rxn_ID>=1))
            assert(all(rxn_ID<=nreactions))
            
        case 'Update_propensities_for_genes'
            gene_IDs=varargin{2};
            
            assert(size(gene_IDs,1)==1) %should be row vector
            assert(size(gene_IDs,2)>=1) %should be at least one gene
            
            assert(all(rem(gene_IDs,1)==0)) %whole numbers
            assert(all(gene_IDs>=1))
            assert(all(gene_IDs<=ngenes))
            
        otherwise
            error('inputs are wrong, see help_state_and_propensities')
    end
            

else
    error('inputs are not what this function expects, see help update_state_and_propensities')
end

%% if update_mode == 'Run_rxn', then change state according to the reaction


if strcmp(update_mode,'Run_rxn')
    
    if rxn_ID<=ngenes %production reaction
        gene_ID_to_update=rxn_ID;
        
        if state(gene_ID_to_update)==1
            error('attempted a production reaction on gene_ID=%i (rxn_ID=%i), but its state was already =1',gene_ID_to_update,rxn_ID)
        end
        
        % update state
        state(gene_ID_to_update)=1;

        
    else % degredation reaction
        gene_ID_to_update=rxn_ID-ngenes; %with 4 genes, rxn_ID=5 means degrade gene 1
        
        if state(gene_ID_to_update)==0
            error('attempted a degredation reaction on gene_ID=%i (rxn_ID=%i), but its state was already =0',gene_ID_to_update,rxn_ID)
        end
        
        %update state
        state(gene_ID_to_update)=0;
        
    end
    
    % now need the gene IDs that need their propensities updated. These
    % include those of the gene that just got updated in
    % update_mode=='Run_rxn' AND the decendents of that gene
    
    % get the successors of the gene whose state we just updated (gene_ID_to_update)
    succInds_gene_to_update=G_Nodes.succInds{gene_ID_to_update}; %row vector, may be empty (1x0 size)
    
    gene_IDs=[gene_ID_to_update,succInds_gene_to_update]; %both the gene and its successors
end

%% Have list of gene_IDs to be updated --> Update it using proMatrix and decMatrix

for indGene_ID=1:length(gene_IDs)
    gene_ID=gene_IDs(indGene_ID);
    
    if state(gene_ID)==0
        %then gene is off--> set decay rate to zero
        propensities(gene_ID+ngenes)=0;
        
        % now figure out what production rate should be. Use proMatrix for
        % this, which is a lookup table for the production rate. The first dimension (Ie. row)
        % is for the first 'proPred' (production predecessor) --> if it is
        % OFF then look in the first row, if it is ON then look in the
        % second row. Similarly, the second dimension (Ie. column) is for
        % the second 'proPred' --> if it is OFF then look in the first column,
        % if it is ON then look in the second column. With more proPred genes, this continues for
        % more dimensions. 
        
        % get proPreds
        proPredInds=G_Nodes.proPredInds{gene_ID}; % [1 x nProPreds]
        proMatrix=G_Nodes.proMatrix{gene_ID};
        nProPreds=size(proPredInds,2);
        
        if nProPreds==0
            propensities(gene_ID)=proMatrix(1); %only one value
        elseif nProPreds==1 %need to get state of production-affecting predecessor
            proPredStates=state(1,proPredInds);
            propensities(gene_ID)=proMatrix(proPredStates+1,1);
        elseif nProPreds==2
            proPredStates=state(1,proPredInds);
            propensities(gene_ID)=proMatrix(proPredStates(1)+1,proPredStates(2)+1);
        else
            error('shouldnt be here, code only meant for 2 proPreds')
        end
    
        
    elseif state(gene_ID)==1
        %then gene is off--> set decay rate to zero
        propensities(gene_ID)=0;
        
        % now figure out what decay rate should be. Use decMatrix. See
        % 'proMatrix' discussion above for details.
        
        decPredInds=G_Nodes.decPredInds{gene_ID}; % [1 x nProPreds]
        decMatrix=G_Nodes.decMatrix{gene_ID};
        nDecPreds=size(decPredInds,2);
        
        if nDecPreds==0
            propensities(gene_ID+ngenes)=decMatrix(1); %only one value
        elseif nDecPreds==1 %need to get state of production-affecting predecessor
            decPredStates=state(1,decPredInds);
            propensities(gene_ID+ngenes)=decMatrix(decPredStates+1,1);
        elseif nDecPreds==2
            decPredStates=state(1,decPredInds);
            propensities(gene_ID+ngenes)=decMatrix(decPredStates(1)+1,decPredStates(2)+1);
        else
            error('shouldnt be here, code only meant for 2 decPreds')
        end
        
    else error('shouldnt be here. state should be 0 or 1')
    end
    
end

state_out=state;
propensities_out=propensities;


end

                                                                                                 %