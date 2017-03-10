for jj = 1:26
    writetable(fitgravcoeff{1,jj+1},'smallgrav3fit_zerorests.xlsx','Sheet',jj,'Range','A1:D5')
end

alldestys = cell(1,11);
alldestid = [];
for ii = 1:size(sortuniqorigins,1)
    % recover the chiefdom id for the origin
    oid = sortuniqorigins(ii);
    opopsy(ii) = pop_counts(oid);
    destys = destinationid(find(oid==originid));
    alldestys{ii} = destys;
    alldestid = [alldestid; destys];
end

%%
fpEb2 = zeros(1,size(windomain,2));
fc2Eb3 = zeros(1,size(windomain,2));
fc2Eb2 = zeros(1,size(windomain,2));
fc3Eb2 = zeros(1,size(windomain,2));
fc3Eb3 = zeros(1,size(windomain,2));
fc3Eb4 = zeros(1,size(windomain,2));

pcutoff = 1e-0;

for ii = 2:size(windomain,2);
    if fitpowercoef{1,ii}.pValue(2)<pcutoff
        fpEb2(ii) = fitpowercoef{1,ii}.Estimate(2);
    end
    if fitgravcoef2{1,ii}.pValue(3)<pcutoff
        fc2Eb3(ii) = fitgravcoef2{1,ii}.Estimate(3);
    end
    if fitgravcoef2{1,ii}.pValue(2)<pcutoff
        fc2Eb2(ii) = fitgravcoef2{1,ii}.Estimate(2);
    end
    if fitgravcoeff{1,ii}.pValue(2)<pcutoff
        fc3Eb2(ii) = fitgravcoeff{1,ii}.Estimate(2);
    end
    if fitgravcoeff{1,ii}.pValue(3)<pcutoff
        fc3Eb3(ii) = fitgravcoeff{1,ii}.Estimate(3);
    end
    if fitgravcoeff{1,ii}.pValue(4)<pcutoff
        fc3Eb4(ii) = fitgravcoeff{1,ii}.Estimate(4);
    end
end
figure; subplot(3,1,1);
plot(windomain,fpEb2); xlabel('time (days after)'); title('fitnlm \tau_1\equiv 0,\tau_2\equiv 0'); ylabel(' \alpha');
subplot(3,1,2);
plot(windomain,fc2Eb3,'DisplayName','\tau_2'); xlabel('time (days after)'); ylabel('\tau'); title('fitnlm \rho\def 0');
hold on; plot(windomain,fc2Eb2,'DisplayName','\tau_1'); xlabel('time (days after)'); title('fitnlm \rho\def 0'); legend toggle
subplot(3,1,3);
plot(windomain,fc3Eb2,'DisplayName','\tau_1'); xlabel('time (days after)'); ylabel('parameter fit'); title('fitnlm');
hold on; 
plot(windomain,fc3Eb3,'DisplayName','\tau_2'); 
plot(windomain,-fc3Eb4,'DisplayName','\rho'); legend toggle

%%
[sortuniqorigins,isort,iorig] = unique(originL3.median);

rests = zeros(size(sortuniqorigins));
moves = zeros(size(sortuniqorigins));
for ii = 1:size(sortuniqorigins,1)
    % recover the chiefdom id for the origin
    oid = sortuniqorigins(ii);
    destys = destinationL3.median(find(oid==originL3.median));
    for jj = 1:size(destys)
        if oid == destys(jj)
            rests(ii) = rests(ii) + 1;
        else
            moves(ii) = moves(ii) + 1;
        end
    end
end

[loma moma] = sort(pop_counts);

touched = find(mean(allrestprob(moma,:),2)>0);
figure; plot(pop_counts(moma(touched)),mean(allrestprob(moma(touched),:),2),'r.')
xlabel('population');
ylabel('mean probability of self-infection');
title('chiefdoms of maximum population');
%%

figure; hold on;
for jj = 1:13
    plot(windomain, allrestprob(moma(touched(jj)),:),'DisplayName', ebolaAdminL3.textAttribute{moma(touched(jj))})
end


%% Drawing color coded maps of probability for spatial models

scalecol = 14;

cmap=jet(180);
cmap = flipud(cmap);

for chid=153
clear gfcolor gmyColor

gfcolor = -log(gravprobnorm(chid,:));
gmyColor = zeros(size(gfcolor,2),3);
gmyColor = cmap(floor(gfcolor*scalecol),:);

ebolaAdminL3=traitVisData('type','AdminL3','textAttribute',regexprep(SierraLeone.AdminL3.NAME_3,' ','') ...
    ,'numberAttribute',newcenters ...
    ,'shapefileID',SierraLeone.AdminL3.ID_3,'cmap',gmyColor,'otherAttribute',newpolygons ...
    ,'otherAttributeType','polygons','parentTextAttribute',SierraLeone.AdminL3.NAME_2 ...
    ,'parentNumberAttribute',false,'parentShapefileID',SierraLeone.AdminL3.ID_2 ...
    ,'parentOtherAttribute',false,'parentOtherAttributeType','AdminL2' ...
    ,'parentType','AdminL2');

figure; title({ebolaAdminL3.textAttribute{chid},'gravity model '});
ebolaAdminL3.plotMap;

clear lfcolor lmyColor

lfcolor = -log(levyprobnorm(chid,:));
lmyColor = zeros(size(lfcolor,2),3);
lmyColor = cmap(floor(lfcolor*scalecol),:);

ebolaAdminL3=traitVisData('type','AdminL3','textAttribute',regexprep(SierraLeone.AdminL3.NAME_3,' ','') ...
    ,'numberAttribute',newcenters ...
    ,'shapefileID',SierraLeone.AdminL3.ID_3,'cmap',lmyColor,'otherAttribute',newpolygons ...
    ,'otherAttributeType','polygons','parentTextAttribute',SierraLeone.AdminL3.NAME_2 ...
    ,'parentNumberAttribute',false,'parentShapefileID',SierraLeone.AdminL3.ID_2 ...
    ,'parentOtherAttribute',false,'parentOtherAttributeType','AdminL2' ...
    ,'parentType','AdminL2');

figure; title({ebolaAdminL3.textAttribute{chid},'Levy flights'});
%ebolaAdminL2.plotMap('colormap',nan(size(ebolaAdminL2.cmap)))
ebolaAdminL3.plotMap;

chid 

end

%%

figure; subplot(2,1,1); hold on;
title('300 - 360 days later median pop, uniform');
xlabel('data point (linkage)');
plot(pgrav,'k--','DisplayName','grav');
plot(plev,'r--','DisplayName','Levy2');
plot(normlev,'r','DisplayName','-ln(N_{Levy2})');
plot(normg,'k','DisplayName','-ln(N_{grav})');
legend toggle

subplot(2,1,2); hold on; 
plot(normg+pgrav,'k','DisplayName', 'loglike grav');
plot(normlev+plev,'r','DisplayName','loglike Levy2');
plot(logLgrav-logLlevy,'b','DisplayName','llgrav-llLevy2');
legend toggle
%%

figure;
title('300 - 360 days later median pop, uniform');
plot(destinationL3_median,logLgrav-logLlevy,'ro','DisplayName','destinations');
ylabel('llratio');

hold on; 
plot(originL3_median,logLgrav-logLlevy,'b+','DisplayName','origins');
ylabel('llratio');
legend toggle

%%

timept = 18;

oL3median = originL3_medianset{timept};
dL3median = destinationL3_medianset{timept};

gd_prob = zeros(size(oL3median));
ld_prob = zeros(size(oL3median));

g_rank = zeros(size(oL3median));
l_rank = zeros(size(oL3median));

for kk=1:size(oL3median,1);
    dL3 = dL3median(kk);
    oL3 = oL3median(kk);
    [sgp, isgp] = sort(gravprobnorm(oL3,:),'descend');
    [lgp, ilgp] = sort(levyprobnorm(oL3,:),'descend');
    g_rank(kk) = find(isgp==dL3);
    l_rank(kk) = find(ilgp==dL3);
    gd_prob(kk) = gravprobnorm(oL3,dL3);
    ld_prob(kk) = levyprobnorm(oL3,dL3);
end

% figure; plot(log(gd_prob./ld_prob));
% figure; hist(log(gd_prob./ld_prob),-2:0.5:4);

figure; 
subplot(3,1,1); hist(log(gd_prob),-10:1:0); title('gravity prob');
subplot(3,1,2); hist(log(ld_prob),-10:1:0); title('levy prob');
subplot(3,1,3); hist(repmat(log(0.5/152),size(oL3median,1),1),-10:1:0); title('uniform');

btuld = sum(ld_prob>0.003289473684211) - sum(ld_prob>0.3);
btugd = sum(gd_prob>0.003289473684211) - sum(gd_prob>0.3);

wtuld = sum(ld_prob<0.003289473684211);
wtugd = sum(gd_prob<0.003289473684211);

btuld/(wtuld+btuld)*100
btugd/(wtugd+btugd)*100


% figure; hist(g_rank-l_rank,[-140:20:20])



