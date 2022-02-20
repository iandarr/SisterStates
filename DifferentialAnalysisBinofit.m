function Meas=DifferentialAnalysisBinofit(T_states,GroupingGenes)

assert(istable(T_states))
%assert(all(ismember(DifferentialAnalysisTypes,{'SisterNaive','SisterControl'})))
DifferentialAnalysisTypes={'SisterNaive','SisterControl'};
assert(all(ismember(GroupingGenes,T_states.Properties.VariableNames)))
if ischar(GroupingGenes)
    GroupingGenes={GroupingGenes};
end

Meas=struct();

allGeneNames=T_states.Properties.VariableNames(~ismember(T_states.Properties.VariableNames,{'cellInd','parentInd'}));
numAllGenes=length(allGeneNames);
numGroupings=length(GroupingGenes);
%T_states

numDifferentialAnalysisTypes=length(DifferentialAnalysisTypes);

for iAnalysis=1:numDifferentialAnalysisTypes
    
    DifferentialAnalysisType=DifferentialAnalysisTypes{iAnalysis};

for iGrouping=1:numGroupings
    groupingGene=GroupingGenes{iGrouping};

    % high for groupingGene
    idxGroupHigh=T_states.(groupingGene)==true;
    % low for groupingGene
    idxGroupLow=T_states.(groupingGene)==false;

    switch DifferentialAnalysisType
        case 'SisterNaive'
            T_states_high=T_states(idxGroupHigh,:);
            T_states_low=T_states(idxGroupLow,:);

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
        otherwise
            error('case not handled')
    end

    % high for groupingGene
    numCellsInHighGroup=height(T_states_high);
    [hiPhatList,hiPciList]=binofit(sum(T_states_high{:,allGeneNames},1),repmat(numCellsInHighGroup,1,numAllGenes));
    Meas.(DifferentialAnalysisType).(groupingGene).hi=[array2table(allGeneNames','VariableNames',{'geneName'}),array2table([hiPhatList',hiPciList,repmat(numCellsInHighGroup,numAllGenes,1)],'VariableNames',{'phat','pciLower','pciUpper','nCells'})];
    % low for groupingGene
    numCellsInLowGroup=height(T_states_low);
    [loPhatList,loPciList]=binofit(sum(T_states_low{:,allGeneNames},1),repmat(numCellsInLowGroup,1,numAllGenes));
    Meas.(DifferentialAnalysisType).(groupingGene).lo=[array2table(allGeneNames','VariableNames',{'geneName'}),array2table([loPhatList',loPciList,repmat(numCellsInLowGroup,numAllGenes,1)],'VariableNames',{'phat','pciLower','pciUpper','nCells'})];

end
end
end
