function G=make_node_rate_matrices(G)
%% compute the  logic matrices that will be used to update reaction rates for each node
% G.Nodes.DecayLogic in based on looking at relationships of input edges
% and the natural decay rate.


% should have a check to make sure these variables are not already in G.Nodes
if  any(ismember(G.Nodes.Properties.VariableNames, {'predNames','predEdgeRows','proPreds','proMatrix','decPreds','decMatrix'}))
    error('make_node_rate_matrices has been called but table G.Nodes already has at least one of the column names that this function will create')
end

nnodes=numnodes(G);
nedges=numedges(G);
G.Nodes.predNames=cell(nnodes,1);
G.Nodes.predInds=cell(nnodes,1);

G.Nodes.predEdgeRows=cell(nnodes,1);
G.Nodes.succNames=cell(nnodes,1);
G.Nodes.succInds=cell(nnodes,1);

for inode=1:nnodes
    NodeName=G.Nodes.Name(inode);
predNames=predecessors(G,NodeName)'; %predecessors, with node names
predInds=predecessors(G,inode)'; %predecessors, with node numbers not names

npredecessors=length(predNames);
succNames=successors(G,NodeName)'; %successors, with node names
succInds=successors(G,inode)'; %successors, with node numbers not names

predEdgeRows=nan(1,npredecessors);
for ipred=1:npredecessors
iPredName=predNames(ipred);
% get the edge 'relationship' property
firstNodeMatches=strcmp(G.Edges.EndNodes(1:nedges,1),iPredName);
lastNodeMatches=strcmp(G.Edges.EndNodes(1:nedges,2),NodeName);
predEdgeRow=find(all([firstNodeMatches,lastNodeMatches],2));
predEdgeRows(1,ipred)=predEdgeRow; 
end

G.Nodes.predNames(inode)={predNames};
G.Nodes.predInds(inode)={predInds};

G.Nodes.predEdgeRows(inode)={predEdgeRows};

G.Nodes.succNames(inode)={succNames};
G.Nodes.succInds(inode)={succInds};
end

%% Now that G.Nodes.predNames and G.Nodes.predEdgeRows are filled in, make a ProMatrix and DecMatrix to figure out the Production and Decay rates based on the precessessor states
% currently only supported for two inputs

G.Nodes.proPreds=cell(nnodes,1);
G.Nodes.proPredInds=cell(nnodes,1);
G.Nodes.proMatrix=cell(nnodes,1);

G.Nodes.decPreds=cell(nnodes,1);
G.Nodes.decPredInds=cell(nnodes,1);
G.Nodes.decMatrix=cell(nnodes,1);

for inode=1:nnodes
    
    predNames=G.Nodes.predNames{inode}; %size 1xN where N is number of predecessors
    predInds=G.Nodes.predInds{inode};
    
    predEdgeRows=G.Nodes.predEdgeRows{inode};
    
    npredecessors=size(predNames,2); % this will be the number of dimensions that the matrix is.
    
    % grab data from the predecessor 'pred' rows about how the relationship
    % affects the production 'pro' or decay 'dec rates. Notably, either
    % through multiplication 'mult' or addition 'add' to the natural rate.
    % if there is multiple edges affecting a node, the order that these are
    % applied to the node will use 'funct_priority'
    rate_affecting=G.Edges.rate_affecting(predEdgeRows);
    funct_type=G.Edges.funct_type(predEdgeRows);
    funct_mag=G.Edges.funct_mag(predEdgeRows);
    funct_priority=G.Edges.funct_priority(predEdgeRows);
    
    %determine whether the predecessors affect production or decay
    isProPred=strcmp(rate_affecting,'pro');
    isDecPred=strcmp(rate_affecting,'dec');

    % check everything either 'pro' or 'dec'
    nProPred=sum(isProPred);
    nDecPred=sum(isDecPred);
    if npredecessors~=nProPred+nDecPred,error('Check that G.Edges.rate_affecting is either pro or deg since npredecessors=%i, nProRelat=%i,nDecRelat=%i',npredecessors,nProPred,nDecPred),end
    indProPred=find(isProPred); %indices of the production predecessors
    indDecPred=find(isDecPred); %indices of the decay predecessors
    
    
    % put these production and decay predecessor names into G.Nodes table
    G.Nodes.proPreds(inode)={predNames(indProPred)};
    G.Nodes.decPreds(inode)={predNames(indDecPred)};
    
    G.Nodes.proPredInds(inode)={predInds(indProPred)};
    G.Nodes.decPredInds(inode)={predInds(indDecPred)};
    
    % now focus on just the production predecessors and get their info
    proFunct_type=funct_type(indProPred); % type is 'mult' or 'add'
    proFunct_mag=funct_mag(indProPred); %magnitude is any number
    proFunct_priority=funct_priority(indProPred);%priority is number or nan. 1 comes first (typically addition) 2 comes second (typically multiplication)
    natProRate=G.Nodes.NatProRate(inode);
    
    % create a matrix that represents what the production rate should be
    % when each of the input variables is a given value
    % production matrix
    if nProPred==0
        G.Nodes.proMatrix(inode)={natProRate}; %singular value, 
    elseif nProPred==1
        
        if strcmp(proFunct_type,'mult')
        temp_proMatrix=natProRate*[1, proFunct_mag]'; %two values
        temp_proMatrix(temp_proMatrix<0)=0;%anything below zero should be zero
        G.Nodes.proMatrix(inode)={temp_proMatrix};
        
        elseif strcmp(proFunct_type,'add')
        temp_proMatrix= natProRate + [0, proFunct_mag]'; %two values
        temp_proMatrix(temp_proMatrix<0)=0;%anything below zero should be zero
        G.Nodes.proMatrix(inode)={temp_proMatrix};
        end

    elseif nProPred>=2   
        error('With current code, cant have >=2 edges that both affect a given rate type (Ie. production or decay). Problem is with node %s',G.Nodes.Name{inode})
        %temp_proMatrix=natProRate*reshape(ones(2^nProPred,1),repmat(2,[1,nProPred]));
        % for later >2 will need to do something like temp_Matrix=natProRate*reshape(ones(2^nProPred,1),repmat(2,[1,nProPred]))
    end
    
    
    
    % now focus on just the decay predecessors and get their info
    decFunct_type=funct_type(indDecPred); % type is 'mult' or 'add'
    decFunct_mag=funct_mag(indDecPred); %magnitude is any number
    decFunct_priority=funct_priority(indDecPred);%priority is number or nan. 1 comes first (typically addition) 2 comes second (typically multiplication)
    natDecRate=G.Nodes.NatDecRate(inode);    
     
     if nDecPred==0
        G.Nodes.decMatrix(inode)={natDecRate}; %singular value, 
    elseif nDecPred==1
        
        if strcmp(decFunct_type,'mult')
        temp_decMatrix=natDecRate*[1, decFunct_mag]'; %two values
        temp_decMatrix(temp_decMatrix<0)=0;%anything below zero should be zero
        G.Nodes.decMatrix(inode)={temp_decMatrix};
        
        elseif strcmp(decFunct_type,'add')
        temp_decMatrix= natDecRate + [0, decFunct_mag]'; %two values
        temp_decMatrix(temp_decMatrix<0)=0;%anything below zero should be zero
        G.Nodes.decMatrix(inode)={temp_decMatrix};
        end

    elseif nDecPred>=2   
        error('With current code, cant have >=2 edges that both affect a given rate type (Ie. production or decay). Problem is with node %s',G.Nodes.Name{inode})
    end
end



end
