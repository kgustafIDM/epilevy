%% hack away all these name differences
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

%%
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

%%
%
% EVDnetworkidSLE(:,1) = nodeL3id{EVDnetworkSLE(:,1)};
% EVDnetworkidSLE(:,2) = nodeL3id{EVDnetworkSLE(:,2)};

%%

daterange = dateboundhigh-dateboundlow;

if strcmp(timertype,'window')
    numtpoints = floor(daterange/twinsize);
    windowcases = cell(1,numtpoints);
    datecutset = 0:twinsize:dateboundhigh-dateboundlow;
elseif strcmp(timertype,'cumulative')
    numtpoints = cumultpts;
    beforecases = cell(1,numtpoints);
    timestep = floor(daterange/numtpoints);
    datecutset = timestep:timestep:timestep*numtpoints;
else
    error('timer type bad');
end

RL1L2set = zeros(1,numtpoints);
Pset = zeros(1,numtpoints);
casesat = zeros(1,numtpoints);
linkagesat = zeros(1,numtpoints);
normRset = zeros(1,numtpoints);
keepidset = cell(1,numtpoints);
originL3_maxset = cell(1,numtpoints);
originL3_minset = cell(1,numtpoints);
originL3_meanset = cell(1,numtpoints);
originL3_medianset = cell(1,numtpoints);
destinationL3_maxset = cell(1,numtpoints);
destinationL3_minset = cell(1,numtpoints);
destinationL3_meanset = cell(1,numtpoints);
destinationL3_medianset = cell(1,numtpoints);

fitgravcoeff = cell(1,numtpoints);
fitgravcoef2 = cell(1,numtpoints);
fitpowercoef = cell(1,numtpoints);
alphaCDF = zeros(1,numtpoints);
likelipl = zeros(1,numtpoints);
coefpowset = zeros(1,numtpoints);

siglinksnetid = find(EVDpvalueSLE<1e-2);

clear datelow datehigh

%%
for k = 1:numtpoints
    k
    datelow = dateboundlow + (k-1)*twinsize;
    datehigh = dateboundlow + k*twinsize;
    if strcmp(timertype,'cumulative')
        datebound = dateboundlow + (k)*timestep;
    end
    
    if strcmp(timertype,'window')
        windowcases{k} = find(ebolaTime<datehigh & ebolaTime>=datelow);
        casesat(k) = size(windowcases{k},1);
        postcasecheck = ismember(EVDlinkagesSLE(:,2),windowcases{k});
    elseif strcmp(timertype,'cumulative')
        beforecases{k} = find(ebolaTime<datebound); % option for cumulative counting
        casesat(k) = size(beforecases{k},1);
        postcasecheck = ismember(EVDlinkagesSLE(:,2),beforecases{k});
    else
        error('bad timertype');
    end
    
    postcaseid = find(postcasecheck>0);
    keepid = intersect(postcaseid,siglinksnetid);
    
    EVDlinkagesSLEkeep = EVDlinkagesSLE(keepid,:);
    linkagesat(k) = size(EVDlinkagesSLEkeep,1);
    %
    % ran gmap_dist.m and saved: save drive_durdist_fullto_kk152.mat drive_dur_seconds_full drive_dist_meters_full
    % load C:\Users\kgustafson\Cell_ParkBedford2015Supp2\seqTrack\drive_durdist_fullto_kk152.mat drive_dur_seconds_full drive_dist_meters_full
    
    EVDroaddist.mean = zeros(size(EVDlinkagesSLEkeep,1),1);
    EVDroaddist.range = zeros(size(EVDlinkagesSLEkeep,1),1);
    EVDroaddist.std = zeros(size(EVDlinkagesSLEkeep,1),1);
    
    EVDroadtime.mean = zeros(size(EVDlinkagesSLEkeep,1),1);
    EVDroadtime.range = zeros(size(EVDlinkagesSLEkeep,1),1);
    EVDroadtime.std = zeros(size(EVDlinkagesSLEkeep,1),1);
    
    EVDpop.mean = zeros(size(EVDlinkagesSLEkeep,1),2);
    EVDpop.range = zeros(size(EVDlinkagesSLEkeep,1),2);
    EVDpop.std = zeros(size(EVDlinkagesSLEkeep,1),2);
    EVDpop.max = zeros(size(EVDlinkagesSLEkeep,1),2);
    EVDpop.min = zeros(size(EVDlinkagesSLEkeep,1),2);
    EVDpop.median = zeros(size(EVDlinkagesSLEkeep,1),2);
    
    originL3_max = zeros(size(EVDlinkagesSLEkeep,1),1);
    originL3_min = zeros(size(EVDlinkagesSLEkeep,1),1);
    originL3_mean = zeros(size(EVDlinkagesSLEkeep,1),1);
    originL3_median = zeros(size(EVDlinkagesSLEkeep,1),1);
    
    destinationL3_max = zeros(size(EVDlinkagesSLEkeep,1),1);
    destinationL3_min = zeros(size(EVDlinkagesSLEkeep,1),1);
    destinationL3_mean = zeros(size(EVDlinkagesSLEkeep,1),1);
    destinationL3_median = zeros(size(EVDlinkagesSLEkeep,1),1);
    
    for ii = 1:size(EVDlinkagesSLEkeep,1)
        nL31 = nodeL3id{EVDlinkagesSLEkeep(ii,1)};
        nL32 = nodeL3id{EVDlinkagesSLEkeep(ii,2)};
        drives = zeros(size(nL31,1)*size(nL32,1),1);
        dtimes = zeros(size(nL31,1)*size(nL32,1),1);
        popnode = zeros(size(nL31,1)*size(nL32,1),2);
        kk = 0;
        for jj = 1:size(nL31,1)
            for ll = 1:size(nL32,1)
                kk = kk + 1;
                drives(kk) = drive_dist_meters_fullsym(nL31(jj),nL32(ll));
                dtimes(kk) = drive_dur_seconds_fullsym(nL31(jj),nL32(ll));
                popnode(kk,1) = SLPpop(find(strcmp(nameL3SLE(nL31(jj)),sortednameL3SLE)));
                popnode(kk,2) = SLPpop(find(strcmp(nameL3SLE(nL32(ll)),sortednameL3SLE)));
            end
        end
        
        EVDroaddist.mean(ii) = mean(drives);
        EVDroadtime.mean(ii) = mean(dtimes);
        EVDpop.mean(ii,:) = mean(popnode,1);
        EVDpop.median(ii,:) = median(popnode,1);
        EVDroaddist.range(ii) = max(drives) - min(drives);
        EVDroadtime.range(ii) = max(dtimes) - min(dtimes);
        EVDpop.range(ii,:) = max(popnode,[],1) - min(popnode,[],1);
        EVDroaddist.std(ii) = std(drives);
        EVDroadtime.std(ii) = std(dtimes);
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
        % THIS is clearly redundant
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
            originL3_max(ii) = nL31;
            originL3_min(ii) = nL31;
            originL3_mean(ii) = nL31;
            originL3_median(ii) = nL31;
        else
            originL3_max(ii) = nL31(omaxchiefdom);
            originL3_min(ii) = nL31(ominchiefdom);
            originL3_mean(ii) = nL31(omeanchiefdom);
            originL3_median(ii) = nL31(omedianchiefdom);
        end
        
        if size(nL32,1) == 1
            destinationL3_max(ii) = nL32;
            destinationL3_min(ii) = nL32;
            destinationL3_mean(ii) = nL32;
            destinationL3_median(ii) = nL32;
        else
            destinationL3_max(ii) = nL32(dmaxchiefdom);
            destinationL3_min(ii) = nL32(dminchiefdom);
            destinationL3_mean(ii) = nL32(dmeanchiefdom);
            destinationL3_median(ii) = nL32(dmedianchiefdom);
        end
        
        %    EVDnetworkdist(ii) = pdist2(latlongL3SLE(ids(1),:),latlongL3SLE(ids(2),:));
    end
    
    % if there are cases in the window or cumulatively
    if casesat(k)>0
        % compute the fits for gravity generalization and Clauset
        fit_gravlevy
        fitgravcoeff{k} = mdlg.Coefficients;
        fitgravcoef2{k} = mdlg2.Coefficients;
        fitpowercoef{k} = mdlpow.Coefficients;
        alphaCDF(k) = aX3cdf;
        likelipl(k) = likliplcdf;
        coefpowset(k) = coefpow(1);
    end
    
    % this function computes the likelihood ratio
    LP_gravlevy
    
    RL1L2set(k) = RL1L2;
    normRset(k) = normRgravlevy;
    Pset(k) = -log10(pgravlevy);
    keepidset{k} = keepid;
    originL3_maxset{k} = originL3_max;
    originL3_minset{k} = originL3_min;
    originL3_meanset{k} = originL3_mean;
    originL3_medianset{k} = originL3_median;
    destinationL3_maxset{k} = destinationL3_max;
    destinationL3_minset{k} = destinationL3_min;
    destinationL3_meanset{k} = destinationL3_mean;
    destinationL3_medianset{k} = destinationL3_median;
    
end
%%
startplotid = find(isnan(normRset(1:end-1)));

normRset(startplotid) = 0;
RL1L2set(startplotid) = 0;
Pset(startplotid) = 0;

figure; subplot(2,1,1); yyaxis left

dName = [dataset,' ',chiefpop];
if strcmp(timertype,'window')
    plot(datecutset(1:end-1)-min(datecutset)+twinsize/2,normRset,'DisplayName',dName); ylabel(['normalized LR ',gravity_type,' ',levyfit]); grid on
elseif strcmp(timertype,'cumulative')
    plot(datecutset,normRset,'DisplayName',dName); ylabel(['normalized LR ',gravity_type,' ',levyfit]); grid on
end
yyaxis right
if strcmp(timertype,'window')
    plot(datecutset(1:end-1)-min(datecutset)+twinsize/2,casesat,'DisplayName',dName); ylabel('case count'); grid on
elseif strcmp(timertype,'cumulative')
    plot(datecutset,casesat,'DisplayName',dName); ylabel('case count'); grid on
end
title(['sequences ',phylotype, ' of linkages N=',num2str(size(EVDlinkagesSLE,1))]);
if strcmp(timertype,'window')
    xlabel(['center of ',num2str(twinsize),'-day window (days after 18-Mar-2014)']);
elseif strcmp(timertype,'cumulative')
    xlabel('days after 18-Mar-2014');
else
    error('timertype wrong');
end
legend toggle;
subplot(2,1,2);
if strcmp(timertype,'window')
    plot(datecutset(1:end-1)-min(datecutset)+twinsize/2,Pset,'DisplayName',dName); ylabel('P score'); grid on
elseif strcmp(timertype,'cumulative')
    plot(datecutset,Pset,'DisplayName',dName); ylabel('P score'); grid on
end
legend toggle;
