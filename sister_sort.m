function S=sister_sort(S,Tin_pop_states,tnow,SortType,twait,sortGeneNames)
% Sout=sister_sort(Sin,Tin_pop_states,'FlagHiLo_Wait_AbortIfChanged'

assert(isstruct(S))
assert(istable(Tin_pop_states))
assert(ismatrix(tnow))
assert(ischar(SortType))
assert(ismatrix(twait))

if strcmp(SortType,'FlagDiff_Wait')
    % flag once sisters have different states of gene . Then, wait twait and
    % sort. Even if sisters are same in the intervening timeperiod, still
    % sort (Ie. don't abort any flagged sisters)
    
elseif strcmp(SortType,'FlagDiff_Wait_AbortIfChanged')
    % flag once sisters have different states of gene __. Then, wait twait and
    % sort. But, if sisters change state in any way during twait, then
    % abort that sort. (Ie. remove flag)
    
else
    error('code not set up to handle other sort types')
end

ncells=height(Tin_pop_states);
isSis1=[true([ncells/2,1]);false([ncells/2,1])];
isSis2=[false([ncells/2,1]);true([ncells/2,1])];
assert(all(Tin_pop_states.parentInd(isSis1)==Tin_pop_states.parentInd(isSis2),1)) %expect top half to be sister1, bottom half sister 2, all in order

assert(all(ismember(sortGeneNames,fieldnames(S)'))) %S must contain a fieldname (referencing a table) for each desired sortGeneName
nSortGenes=size(sortGeneNames,1);

% get all gene names which, will will be needed when we sort cells. Base
% this off table VariableNames
isGeneNameInd=~ismember(Tin_pop_states.Properties.VariableNames,{'cellInd','parentInd'});
allGeneNames=Tin_pop_states.Properties.VariableNames(isGeneNameInd);

for indGene=1:nSortGenes
    % check that Sin has all the appropriate fieldnames
    geneName=sortGeneNames{indGene};
    if ~all(ismember({'cellInd','parentInd','FlagTime','FlagType','SortTime'},S.(geneName).Properties.VariableNames')),
        error('asked to sort but at least one table in structure Sin doesnt have required fieldnames')
    elseif ~all(ismember(geneName,S.(geneName).Properties.VariableNames'))
        error('the gene that were sorting off of is not a VariableName in at least one table in the Sin structure')
    end


end


for indGene=1:nSortGenes
    geneName=sortGeneNames{indGene};
    
    isFlagged=~isnan(S.(geneName).FlagTime);
    isSorted=~isnan(S.(geneName).SortTime);
    
    isFlagTypeHi=strcmp(S.(geneName).FlagType,'hi');
    isFlagTypeLo=strcmp(S.(geneName).FlagType,'lo');
    
    assert(all(isFlagged==sum([isFlagTypeHi,isFlagTypeLo],2)))

    %% Find any sisters to abort a previous flagging (Ie. remove the flag). If this is turned on.
    if strcmp(SortType,'FlagDiff_Wait_AbortIfChanged')
    
    % check if the gene states are still what is expected
    isChangedWasFlaggedHi=all([isFlagTypeHi,~isSorted,1~=Tin_pop_states.(geneName)],2);
    isChangedWasFlaggedLo=all([isFlagTypeLo,~isSorted,0~=Tin_pop_states.(geneName)],2);
    
    isChangedWasFlagged=any([isChangedWasFlaggedHi,isChangedWasFlaggedLo],2);
    % but this may only have been 1 cell out of the pair that changed. If
    % it was indeed only one of the sisters that changed, identify that
    % sister too:
    isSisOfChanged=[isChangedWasFlagged(isSis2);isChangedWasFlagged(isSis1)];
    
    isChangedOrSisOfChanged=any([isChangedWasFlagged,isSisOfChanged],2);
    
    
    % remove flag from these changed cells
    S.(geneName).FlagTime(isChangedOrSisOfChanged)=nan;
    S.(geneName).FlagType(isChangedOrSisOfChanged)=cell(sum(isChangedOrSisOfChanged),1);
    end
    
    %% Next, find any sister cells that should be newly flagged (Ie. add 'FlagTime' and 'FlagType' in S)
    sisIsDiff=Tin_pop_states.(geneName)(isSis1)~=Tin_pop_states.(geneName)(isSis2); %sister is different
    isDiff=[sisIsDiff;sisIsDiff];% back to whole population indices;
    
    state=Tin_pop_states.(geneName);
    
    isToBeFlaggedHi=all([~isSorted,~isFlagged,isDiff,state==1],2);
    isToBeFlaggedLo=all([~isSorted,~isFlagged,isDiff,state==0],2);
    
    % now updated FlagTime and SortType for these
    S.(geneName).FlagTime(any([isToBeFlaggedHi,isToBeFlaggedLo],2))=tnow;
    S.(geneName).FlagType(isToBeFlaggedHi)={'hi'};
    S.(geneName).FlagType(isToBeFlaggedLo)={'lo'};
    
    %% Next, sort out any cells where we've waited long enough after the flag (Ie. update S with 'SortTime' and record all gene states)
    
    hasWaitedLongEnough=tnow - S.(geneName).FlagTime>=twait;
    isToBeSorted=all([hasWaitedLongEnough,~isSorted],2);
    
    S.(geneName).SortTime(isToBeSorted)=tnow;
    S.(geneName)(isToBeSorted,allGeneNames)=Tin_pop_states(isToBeSorted,allGeneNames);
    
end


end