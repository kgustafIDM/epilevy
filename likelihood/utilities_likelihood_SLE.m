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

for ii = 5:size(windomain,2);
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
%%
cmap = jet(15);
cmapg = gray(15);

[ebTsort, st2] = sort(ebolaTree.time);

origindates = ebolaTree.time(ebolaTree.network(:,1));
destinationdates = ebolaTree.time(ebolaTree.network(:,2));

[sortodates, soids] = sort(origindates);
uniqodates = unique(sortodates);

nodeL2sort = nodeL2id(st2);

mintime = 735676;
figure; hold on;
for kk = 1:size(uniqodates,1)
    currentid = find(origindates==uniqodates(kk));
    for jj=1:size(currentid,1)
        XX = [origindates(currentid(jj))-mintime,destinationdates(currentid(jj))-mintime];
        YY = [nodeL2id(ebolaTree.network(currentid(jj),1)),nodeL2id(ebolaTree.network(currentid(jj),2))];
        plot(XX,YY,'Color',cmapg(10,:));
        scatter(XX(1),YY(1),'MarkerEdgeColor',cmap(YY(1),:));
        scatter(XX(2),YY(2),'.','MarkerFaceColor',cmap(YY(2),:));
        ax = gca; ax.YTick = (1:14); ax.YTickLabel = nameL2SLE;
    end
    pause(0.2); drawnow;
end

%%
cmap=bone;
chid = 97;

f = gravprobnorm(chid,:);
f = f*floor(1/max(f));

cm = cmap; % returns the current color map
myColor = zeros(size(f,2),3);
for jj = 1:size(f,2)
    colorID = max(1, sum(f(jj) > [0:1/length(cm(:,1)):1]));
    myColor(jj,:) = cm(colorID, :); % returns your color
end

ebolaAdminL3=traitVisData('type','AdminL3','textAttribute',regexprep(SierraLeone.AdminL3.NAME_3,' ',''),'numberAttribute',newcenters ...
    ,'shapefileID',SierraLeone.AdminL3.ID_3,'cmap',myColor,'otherAttribute',newpolygons ...
    ,'otherAttributeType','polygons','parentTextAttribute',SierraLeone.AdminL3.NAME_2 ...
    ,'parentNumberAttribute',false,'parentShapefileID',SierraLeone.AdminL3.ID_2 ...
    ,'parentOtherAttribute',false,'parentOtherAttributeType','AdminL2' ...
    ,'parentType','AdminL2');
%ebolaAdminL3.createColormap('Kakua','Nigeria');

% [ix,loc]=ismember({'kissiteng','kissitongi','jawie','mambolo','kissikama','luawa','njaluahun','pehebongre','malema','kakua','pejewest','nongowa'},lower(ebolaAdminL3.textAttribute));
% cmap2=cbrewer('qual','Set1',length(loc));
% cmap2(3,:)=[0, 1,0];
% cmap2(8,:)=[0,.9,.9];
% cmap=ebolaAdminL3.cmap;
% sortIdx=[1,3,4,5,6,7,8,2,9,10,11,12];
% cmap(loc(sortIdx),:)=cmap2;
% ebolaAdminL3.createColormap(cmap);

% figure
% ebolaAdminL2.plotMap('colormap',nan(size(ebolaAdminL2.cmap)))
% ebolaAdminL3.plotMap('includeonly',{'kissiteng','jawie','luawa','mambolo','freetown1','valunia'})
figure; subplot(2,1,1); title(['gravity model ',ebolaAdminL3.textAttribute{chid}]);
%ebolaAdminL2.plotMap('colormap',nan(size(ebolaAdminL2.cmap)))
ebolaAdminL3.plotMap;

f = levyprobnorm(chid,:);
f = f*floor(1/max(f));

cm = cmap; % returns the current color map
myColor = zeros(size(f,2),3);
for jj = 1:size(f,2)
    colorID = max(1, sum(f(jj) > [0:1/length(cm(:,1)):1]));
    myColor(jj,:) = cm(colorID, :); % returns your color
end

ebolaAdminL3.cmap = myColor;

subplot(2,1,2); title('Levy flights');
%ebolaAdminL2.plotMap('colormap',nan(size(ebolaAdminL2.cmap)))
ebolaAdminL3.plotMap;


%%

figure; subplot(2,1,1); hold on;
title('330 - 390 days later median pop, uniform');
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
