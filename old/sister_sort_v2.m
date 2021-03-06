function S=sister_sort_v2(S,Tin_pop_states,tnow) %,SortType,twait,sortGeneNames)
% Sout=sister_sort(Sin,Tin_pop_states,'FlagHiLo_Wait_AbortIfChanged'
% v2 want to include support for:
%   S to be a cell array of structures
%   S has a field settings as inputs
%   now can specify up to 3 flags
%   now it expects Tin_pop_states to have sisters in adjacent rows (Ie.
%   rows 1 and 2 are sisters, 3 and 4 are sisters, etc...)

% S(1).settings.sortGeneNames=G.Nodes.Name;
% S(1).settings.Flag1='lo-lo';
% S(1).settings.Flag2='hi-lo';
% S(1).settings.Flag3=[];
% S(1).settings.SortTrigger='Flag2';
% S(1).settings.WaitTimeAfterSortTrigger=1;
% S(1).settings.AbortEvent=[];
% assert(isstruct(S))
% assert(istable(Tin_pop_states))
% assert(ismatrix(tnow))
% assert(ischar(SortType))
% assert(ismatrix(twait))

%% check over Tin_pop_states and check over 

ncells=height(Tin_pop_states);
assert(rem(ncells,2)==0) %must be even number

isSis1=repmat([true;false],[ncells/2 1]);
isSis2=repmat([false;true],[ncells/2 1]);
assert(all(Tin_pop_states.parentInd(isSis1)==Tin_pop_states.parentInd(isSis2),1)) %expect top half to be sister1, bottom half sister 2, all in order

indOfSisterLookup=[[1:ncells]' + repmat([1 -1]',[ncells/2 1])];
%% Now look at sort structure S.

if ~isstruct(S),error('S should be a structure'),end

nSortStructs=length(S);

if ~any(strcmp(fieldnames(S),'settings'))
    error('you must input S with a settings fieldname')
end

if ~ismember({'data'},fieldnames(S))
    S(1).data=[];
end


for iS=1:nSortStructs
    %% Look over expected input fields to S
    
    expected_fields={'sortGeneNames','Flag1','Flag2','Flag3','SortTrigger','WaitTimeAfterSortTrigger','AbortEvent'};
    missing_field_ind=find(~ismember(expected_fields,fieldnames(S(iS).settings)));
    if ~isempty(missing_field_ind)
        error(['these field(s) are missing from S(%i).settings:\n',repmat('%s\n',[1,length(missing_field_ind)])],iS,expected_fields{missing_field_ind})
    end
    
    sortGeneNames=S(iS).settings.sortGeneNames;
    
    isGeneNameInd=~ismember(Tin_pop_states.Properties.VariableNames,{'cellInd','parentInd'});
    allGeneNames=Tin_pop_states.Properties.VariableNames(isGeneNameInd);
    
    %% if called for first time, need to fill in S(iS).data
    %then check if its empty
    if isempty(S(iS).data) % then still initialize S(iS).data
        flag_ID_values=[NaN false true];
        flag_ID_categories={'-','lo','hi'};
     
        for igene=1:length(sortGeneNames)
            sortGeneName=sortGeneNames{igene};
            S(iS).data.(sortGeneName)=[...
                Tin_pop_states(:,{'cellInd','parentInd'}),...
                array2table(nan(ncells,1),'VariableNames',{'Flag1_time'}),...
                array2table(categorical(nan(ncells,1),flag_ID_values,flag_ID_categories),'VariableNames',{'Flag1_ID'}),...
                array2table(nan(ncells,1),'VariableNames',{'Flag2_time'}),...
                array2table(categorical(nan(ncells,1),flag_ID_values,flag_ID_categories),'VariableNames',{'Flag2_ID'}),...
                array2table(nan(ncells,1),'VariableNames',{'Flag3_time'}),...
                array2table(categorical(nan(ncells,1),flag_ID_values,flag_ID_categories),'VariableNames',{'Flag3_ID'}),...
                array2table(nan(ncells,1),'VariableNames',{'Sort_time'}),...
                array2table(nan(ncells,1),'VariableNames',{'sortGroup'}),...
                array2table(false(ncells,length(sortGeneNames)),'VariableNames',sortGeneNames)];
        end
        
        
    end %end initialization of sort tables within S(iS).data
    
    %% now just make sure the genes that we'll be using to sort based on are
    % actually named in the structure.
    assert(all(ismember(sortGeneNames,fieldnames(S(iS).data)))) %S.data (a table) must contain a VariableName for each desired sortGeneName
    
    nSortGenes=length(sortGeneNames);
    
    % this part hast to be done as of 28apr19
    for iSortGene=1:nSortGenes
        
        %% get general information about the states
        sortGeneName=sortGeneNames{iSortGene};
        twait=S(iS).settings.WaitTimeAfterSortTrigger;
        
        isHi=S(iS).data.(sortGeneName).(sortGeneName)==true;
        isLo=S(iS).data.(sortGeneName).(sortGeneName)==false;
        
        bothSisAreLo=all([isLo(isSis1),isLo(isSis2)],2); % only is ncell/2 height
        bothSisAreLo=reshape([bothSisAreLo';bothSisAreLo'],[length(bothSisAreLo)*2,1]);
  
        bothSisAreHi=all([isHi(isSis1),isHi(isSis2)],2); %only ncell/2 height
        bothSisAreHi=reshape([bothSisAreHi';bothSisAreHi'],[length(bothSisAreHi)*2,1]);
        
        sistersAreDiffStates=all([~bothSisAreLo,~bothSisAreHi],2);
        
        %assert(all(sum([isHi,isLo],2)==1))
        isSorted=~isnan(S(iS).data.(sortGeneName).Sort_time);
        %% Unflagged --> Flag1
        % Apply Flag1 to all unflagged that should get it
        isEligibleForFlag1=~isnan(S(iS).data.(sortGeneName).Flag1_time);
        assert(all(isnan(S(iS).data.(sortGeneName).Sort_time(isEligibleForFlag1))))
        assert(all(isnan(S(iS).data.(sortGeneName).Flag2_time(isEligibleForFlag1))))
        if strcmp(S(iS).settings.Flag1,'lo-lo')  
            makeFlag1=all([isEligibleForFlag1,bothSisAreLo],2);
        elseif strcmp(S(iS).settings.Flag1,'hi-hi')
            makeFlag1=all([isEligibleForFlag1,bothSisAreHi],2);
        else
            error('didnt think of this')
        end
        
        makeFlag1_IDLo=all([makeFlag1,isLo],2);
        makeFlag1_IDHi=all([makeFlag1,isHi],2);
        S(iS).data.(sortGeneName).Flag1_time(makeFlag1)=tnow;
        S(iS).data.(sortGeneName).Flag1_ID(makeFlag1_IDLo)='lo';
        S(iS).data.(sortGeneName).Flag1_ID(makeFlag1_IDHi)='hi';
        
        %% Remove Flag1 (will not occur at this point, can skip this)
        
        %% Flag1 --> Flag2
        % Next, find any sister cells that should be newly flagged (Ie. add 'FlagTime' and 'FlagType' in S)
        isEligibleForFlag2=~isnan(S(iS).data.(sortGeneName).Flag2_time);
        assert(all(isnan(S(iS).data.(sortGeneName).Sort_time(isEligibleForFlag2))))
        
        if ismember(S(iS).settings.Flag2,{'lo-hi','hi-lo'})
            makeFlag2=all([isEligibleForFlag2,sistersAreDiffStates],2);
        else
            error('didnt think of this')
        end

        makeFlag2_IDLo=all([makeFlag2,isLo],2);
        makeFlag2_IDHi=all([makeFlag2,isHi],2);
        S(iS).data.(sortGeneName).Flag2_time(makeFlag2)=tnow;
        S(iS).data.(sortGeneName).Flag2_ID(makeFlag2_IDLo)='lo';
        S(iS).data.(sortGeneName).Flag2_ID(makeFlag2_IDHi)='hi';

        %% Flag2 --> Unflagged, if cells are different now
        if strcmp(S(iS).settings.AbortEvent,'AbortSortIfChangedWhileWaiting')
            isFlag2=all([~isSorted,~isnan(S(iS).data.(sortGeneName).Flag2_time)],2);
            isFlag2hi=all([isFlag2,S(iS).data.(sortGeneName).Flag2_ID=='hi'],2);
            isFlag2lo=all([isFlag2,S(iS).data.(sortGeneName).Flag2_ID=='lo'],2);
            flag2HasChanged=any([all([isFlag2hi,~isHi],2),all([isFlag2lo,~isLo],2)],2);
            indFlag2HasChanged=find(flag2HasChanged);
            indRemoveFlag2=[indFlag2HasChanged;indOfSisterLookup(flag2HasChanged)];
            S(iS).data.(sortGeneName).Flag2_time(indRemoveFlag2)=nan;
            S(iS).data.(sortGeneName).Flag2_ID(indRemoveFlag2)='-';
            S(iS).data.(sortGeneName).Flag1_time(indRemoveFlag2)=nan;
            S(iS).data.(sortGeneName).Flag1_ID(indRemoveFlag2)='-';
        end 
        
        %% Flag2 --> Sort (if twait has passed)
        % Next, sort out any cells where we've waited long enough after the flag (Ie. update S with 'SortTime' and record all gene states)

        isToBeSorted=all([~isSorted,tnow - S(iS).data.(sortGeneName).Flag2_time>=twait],2);
        S(iS).data.(sortGeneName).Sort_time(isToBeSorted)=tnow;
        S(iS).data.(sortGeneName).sortGroup(all([isToBeSorted,isHi],2))=true;
        S(iS).data.(sortGeneName).sortGroup(all([isToBeSorted,isLo],2))=false;
        S(iS).data.(sortGeneName)(isToBeSorted,allGeneNames)=Tin_pop_states(isToBeSorted,allGeneNames);

    end % each sortGeneName
    
    
end %each iS

end % end this S(iS)