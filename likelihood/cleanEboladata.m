
if needtoload
    addpath U:\EpiLevy\code_base\epilevy\
    if strcmp(dataset,'Park')
        Park_SLe_network_prelims
    elseif strcmp(dataset,'Rambaut')
        Rambaut_SLe_network_prelims
    else
        error('nodataset');
    end
end

if strcmp(phylotype,'tree')
    nodenameL2 = ebolaTree.sourceData.AdminL2;
    if strcmp(dataset,'Rambaut')
        load nodenameL3_RambautwithPark nodenameL3RP
        nodenameL3 = nodenameL3RP;
    else
        nodenameL3 = ebolaTree.sourceData.AdminL3;
    end
    ebolaTime = ebolaTree.time;
    EVDlinkagesSLE = ebolaTree.network;
    EVDpvalueSLE = ebolaTree.p;
    EVDlink_dur = ebolaTree.duration;
elseif strcmp(phylotype,'network')
    nodenameL2 = ebolaNetwork.sourceData.AdminL2;
    nodenameL3 = ebolaNetwork.sourceData.AdminL3;
    ebolaTime = ebolaNetwork.time;
    EVDlinkagesSLE = ebolaNetwork.network;
    EVDpvalueSLE = ebolaNetwork.p;
    EVDlink_dur = ebolaNetwork.duration;
end

% hack away all these name differences
if strcmp(dataset,'Rambaut')
    idsFreetown = find(strcmp(nodenameL2,'Freetown'));
    for kk = 1:size(idsFreetown,1)
        nodenameL2{idsFreetown(kk)}='WesternUrban';
    end
    idWU = find(strcmp(nodenameL2,'WesternArea'));
    for kk = 1:size(idWU,1)
        nodenameL2{idWU(kk)}='WesternUrban';
    end
elseif strcmp(dataset,'Park')
    nodenameL2{strcmp(nodenameL2,'Koinodugu')}='Koinadugu';
    nodenameL2{strcmp(nodenameL2,'Portloko')}='PortLoko';
    %nodenameL2{strcmp(nodenameL2,'WesternArea')}='WesternUrban';
    idWU = find(strcmp(nodenameL2,'WesternArea'));
    for kk = 1:size(idWU,1)
        nodenameL2{idWU(kk)}='WesternUrban';
    end
    idsFreetown = find(strcmp(nodenameL2,'Freetown'));
    for kk = 1:size(idsFreetown,1)
        nodenameL2{idsFreetown(kk)}='WesternUrban';
    end
else
    error('no dataset');
end

% this is for removing the cases from Guinea, source of importation
nodeL2id = [];
for k = 1:size(nodenameL3,1)
    if sum(find(k==[88 89 90]))==0
        nodeL2id(k) = find(strcmp(nodenameL2(k),nameL2SLE));
    else
        nodeL2id(k) = 1;
    end
end


% here we account for those cases that have only reported AdminL2, while
% keeping the higher accuracy AdminL3 reported cases - an array of chiefdom
% indices is stored for each district, so that the mean road travel distance
% between any of the chiefdoms in the district will be used to represent the
% distance between cases. Maybe this is a bit weird, and one should use the
% district centroid locations, but I had already computed the driving
% distances for each chiefdom to all other chiefdoms

nodeL3id = cell(size(nodenameL3));
for k = 1:size(nodenameL3,1)
    if strcmp(nodenameL3(k),'unknown')==1
        nodeL3id{k} = find(eAdL3inL2==nodeL2id(k));
    elseif strcmp(nodenameL3(k),'Other')==1
        nodeL3id{k} = find(eAdL3inL2==nodeL2id(k));
    else
        nodeL3id{k} = find(strcmp(nodenameL3(k),nameL3SLE));
    end
end
