function Meas=DifferentialAnalysisBinofit2(T_states,idxHiAtSortTime,randSeed)
if nargin>=3
    rng(randSeed)
end
% DifferentialAnalysisBinofit2 now accepts a list of indices
% (idxHiAtSortTime) instead of the sort genes, since with this function version the current state is
% not the state at time of sort

assert(istable(T_states))
DifferentialAnalysisTypes={'SisterControl','SisterNaive','SisterNaiveDownsampled'};


Meas=struct();

allGeneNames=T_states.Properties.VariableNames(~ismember(T_states.Properties.VariableNames,{'cellInd','parentInd'}));
numAllGenes=length(allGeneNames);
%T_states

numDifferentialAnalysisTypes=length(DifferentialAnalysisTypes);

for iAnalysis=1:numDifferentialAnalysisTypes
    
    DifferentialAnalysisType=DifferentialAnalysisTypes{iAnalysis};

    % high for sortGene at sort time
    idxGroupHigh= idxHiAtSortTime;
    
    % low for sortGene at sort time
    idxGroupLow=~idxHiAtSortTime;

    switch DifferentialAnalysisType

        case 'SisterControl'

            cellIndGroupHigh=T_states.cellInd(idxGroupHigh);
            parentIndGroupHigh=T_states.parentInd(idxGroupHigh);
            [~,uniqueParentIndInd,ic]=unique(parentIndGroupHigh);

            indIndParent1HiSibling=uniqueParentIndInd(1==sum(ic==[1:max(ic)],1));
            parentInd1HiSibling=parentIndGroupHigh(indIndParent1HiSibling);

            idxHiAndDifferentFromSibling=idxGroupHigh & any(T_states.parentInd==parentInd1HiSibling',2);
            idxLoAndDifferentFromSibling=idxGroupLow  & any(T_states.parentInd==parentInd1HiSibling',2);

            T_states_high=T_states(idxHiAndDifferentFromSibling,:);
            T_states_low=T_states(idxLoAndDifferentFromSibling,:);
            assert(height(T_states_high)==height(T_states_low))
            numSisterControlCellsPerGroup=height(T_states_high);
            
        case 'SisterNaive'
            T_states_high=T_states(idxGroupHigh,:);
            T_states_low=T_states(idxGroupLow,:);

        case 'SisterNaiveDownsampled'
            indGroupHigh=find(idxGroupHigh);
            indGroupLow=find(idxGroupLow);
            numCellsGroupHigh=sum(idxGroupHigh);
            numCellsGroupLow=sum(idxGroupLow);
            numCellsToSample=min([numCellsGroupHigh,numCellsGroupLow,numSisterControlCellsPerGroup]);
            % I do not set random seed with rng here
            indGroupHighDownsampled=indGroupHigh(randperm(numCellsGroupHigh,numCellsToSample));
            indGroupLowDownsampled=indGroupLow(randperm(numCellsGroupLow,numCellsToSample));
            T_states_high=T_states(indGroupHighDownsampled,:);
            T_states_low=T_states(indGroupLowDownsampled,:);
        otherwise
            error('case not handled')
    end

    % high for sortGene at sort time
    numCellsInHighGroup=height(T_states_high);
    [hiPhatList,hiPciList]=binofit(sum(T_states_high{:,allGeneNames},1),repmat(numCellsInHighGroup,1,numAllGenes));
    Meas.(DifferentialAnalysisType).hi=[array2table(allGeneNames','VariableNames',{'geneName'}),array2table([hiPhatList',hiPciList,repmat(numCellsInHighGroup,numAllGenes,1)],'VariableNames',{'phat','pciLower','pciUpper','nCells'})];
    % low for sortGene at sort time
    numCellsInLowGroup=height(T_states_low);
    [loPhatList,loPciList]=binofit(sum(T_states_low{:,allGeneNames},1),repmat(numCellsInLowGroup,1,numAllGenes));
    Meas.(DifferentialAnalysisType).lo=[array2table(allGeneNames','VariableNames',{'geneName'}),array2table([loPhatList',loPciList,repmat(numCellsInLowGroup,numAllGenes,1)],'VariableNames',{'phat','pciLower','pciUpper','nCells'})];

end
end
