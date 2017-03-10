function [oL3,dL3,EVDpop,EVDroaddist] = get_orig_dest(EVDlink,dd_meters_fullsym,dd_seconds_fullsym,nL3id,SLEpop,nL3SLE,snL3SLE)

EVDroaddist.mean = zeros(size(EVDlink,1),1);
EVDroaddist.median = zeros(size(EVDlink,1),1);
EVDroaddist.max = zeros(size(EVDlink,1),1);
EVDroaddist.min = zeros(size(EVDlink,1),1);
EVDroaddist.range = zeros(size(EVDlink,1),1);
EVDroaddist.std = zeros(size(EVDlink,1),1);

EVDroadtime.mean = zeros(size(EVDlink,1),1);
EVDroadtime.range = zeros(size(EVDlink,1),1);
EVDroadtime.std = zeros(size(EVDlink,1),1);

EVDpop.mean = zeros(size(EVDlink,1),2);
EVDpop.range = zeros(size(EVDlink,1),2);
EVDpop.std = zeros(size(EVDlink,1),2);
EVDpop.max = zeros(size(EVDlink,1),2);
EVDpop.min = zeros(size(EVDlink,1),2);
EVDpop.median = zeros(size(EVDlink,1),2);

oL3.max = zeros(size(EVDlink,1),1);
oL3.min = zeros(size(EVDlink,1),1);
oL3.mean = zeros(size(EVDlink,1),1);
oL3.median = zeros(size(EVDlink,1),1);

dL3.max = zeros(size(EVDlink,1),1);
dL3.min = zeros(size(EVDlink,1),1);
dL3.mean = zeros(size(EVDlink,1),1);
dL3.median = zeros(size(EVDlink,1),1);

for ii = 1:size(EVDlink,1)
    nL31 = nL3id{EVDlink(ii,1)};
    nL32 = nL3id{EVDlink(ii,2)};
    drives = zeros(size(nL31,1)*size(nL32,1),1);
    dtimes = zeros(size(nL31,1)*size(nL32,1),1);
    popnode = zeros(size(nL31,1)*size(nL32,1),2);
    kk = 0;
    for jj = 1:size(nL31,1)
        for ll = 1:size(nL32,1)
            kk = kk + 1;
            drives(kk) = dd_meters_fullsym(nL31(jj),nL32(ll));
            dtimes(kk) = dd_seconds_fullsym(nL31(jj),nL32(ll));
            popnode(kk,1) = SLEpop(find(strcmp(nL3SLE(nL31(jj)),snL3SLE)));
            popnode(kk,2) = SLEpop(find(strcmp(nL3SLE(nL32(ll)),snL3SLE)));
        end
    end
    
    EVDroaddist.mean(ii) = mean(drives);
    EVDroaddist.median(ii) = median(drives);
    EVDroaddist.std(ii) = std(drives);
    EVDroaddist.max(ii) = max(drives);
    EVDroaddist.min(ii) = min(drives);
    EVDroaddist.range(ii) = max(drives) - min(drives);
    
    EVDroadtime.mean(ii) = mean(dtimes);
    EVDroadtime.range(ii) = max(dtimes) - min(dtimes);
    EVDroadtime.std(ii) = std(dtimes);
    
    EVDpop.mean(ii,:) = mean(popnode,1);
    EVDpop.median(ii,:) = median(popnode,1);
    EVDpop.range(ii,:) = max(popnode,[],1) - min(popnode,[],1);
    EVDpop.std(ii,:) = std(popnode,0,1);
    
    % find the chiefdom with the maximum population in the district
    omaxchiefdom = find(max(popnode(:,1))==unique(popnode(:,1),'stable'));
    dmaxchiefdom = find(max(popnode(:,2))==unique(popnode(:,2),'stable'));
    EVDpop.max(ii,1) = max(popnode(:,1));
    EVDpop.max(ii,2) = max(popnode(:,2));
    
    % find the chiefdom with population nearest to the mean population in the district
    val = EVDpop.mean(ii,1); %value to find
    tmp = abs(unique(popnode(:,1),'stable') - val);
    [omed, omeanchiefdom] = min(tmp); %index of closest value
    val = EVDpop.mean(ii,2); %value to find
    tmp = abs(unique(popnode(:,2),'stable') - val);
    [dmed, dmeanchiefdom] = min(tmp); %index of closest value
    
    % find the chiefdom with population nearest to the median population in the district
    % not redundant when there are an even number of chiefdoms in a district
    val = EVDpop.median(ii,1); %value to find
    tmp = abs(unique(popnode(:,1),'stable') - val);
    [omdd, omedianchiefdom] = min(tmp); %index of closest value
    val = EVDpop.median(ii,2); %value to find
    tmp = abs(unique(popnode(:,2),'stable') - val);
    [dmdd, dmedianchiefdom] = min(tmp); %index of closest value
    
    EVDpop.max(ii,1) = max(popnode(:,1));
    EVDpop.max(ii,2) = max(popnode(:,2));
    
    ominchiefdom = find(min(popnode(:,1))==unique(popnode(:,1),'stable'));
    dminchiefdom = find(min(popnode(:,2))==unique(popnode(:,2),'stable'));
    EVDpop.min(ii,1) = min(popnode(:,1));
    EVDpop.min(ii,2) = min(popnode(:,2));
    
    if size(nL31,1) == 1
        oL3.max(ii) = nL31;
        oL3.min(ii) = nL31;
        oL3.mean(ii) = nL31;
        oL3.median(ii) = nL31;
    else
        oL3.max(ii) = nL31(omaxchiefdom);
        oL3.min(ii) = nL31(ominchiefdom);
        oL3.mean(ii) = nL31(omeanchiefdom);
        oL3.median(ii) = nL31(omedianchiefdom);
    end
    
    if size(nL32,1) == 1
        dL3.max(ii) = nL32;
        dL3.min(ii) = nL32;
        dL3.mean(ii) = nL32;
        dL3.median(ii) = nL32;
    else
        dL3.max(ii) = nL32(dmaxchiefdom);
        dL3.min(ii) = nL32(dminchiefdom);
        dL3.mean(ii) = nL32(dmeanchiefdom);
        dL3.median(ii) = nL32(dmedianchiefdom);
    end
    
    %    EVDnetworkdist(ii) = pdist2(latlongL3SLE(ids(1),:),latlongL3SLE(ids(2),:));
end

end