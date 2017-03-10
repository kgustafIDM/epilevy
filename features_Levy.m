
kk = 0;
clearcases = [];    
for ii = 1:6
    for jj = [1:3,9:11]
        kk = kk + 1;
        clearcases(kk) = jj+11*(ii-1);
    end
end

clearsub = [];
clearsup = [];
kk = 0;

for ii = 1:6
    for jj = 1:3
        kk = kk + 1;
        clearsub(kk) = jj+11*(ii-1);
    end
end

kk = 0;
for ii = 1:6
    for jj = 9:11
        kk = kk + 1;
        clearsup(kk) = jj+11*(ii-1);
    end
end

pcc = pdf_slopes(clearcases,:);
dcc = disp_slopes(clearcases,:);
fcc = fp_binary(clearcases,:);

pcsub = pdf_slopes(clearsub,:);
dcsub = disp_slopes(clearsub,:);
fcsub = fp_binary(clearsub,:);

pcsup = pdf_slopes(clearsup,:);
dcsup = disp_slopes(clearsup,:);
fcsup = fp_binary(clearsup,:);

clearsubcases = [1 2 3 7 8 9 13 14 15 19 20 21 25 26 27 31 32 33];
clearsupcases = [4 5 6 10 11 12 16 17 18 22 23 24 28 29 30 34 35 36];

% concatenate features matrix
features_Levyclear = [pcc dcc fcc];
% compute cross-correlation for exploration
CC = corr(features_Levyclear,features_Levyclear);
% compute normalization factor
% mean
mean_featLevy = repmat(mean(features_Levyclear),size(features_Levyclear,1),1);
% variance
var_featLevy = repmat(var(features_Levyclear),size(features_Levyclear,1),1);
norm_featLevy = (features_Levyclear-mean_featLevy)./var_featLevy;

% quick PCA analysis
%[wcoeff,score,latent,tsquared,explained] = pca(features_Levy,'VariableWeights',1./nfact_featLevy);
[wcoeff,score,latent,tsquared,explained] = pca(norm_featLevy);

figure; plot3(score(clearsubcases,1),score(clearsubcases,2),score(clearsubcases,3),'+');
hold on; plot3(score(clearsupcases,1),score(clearsupcases,2),score(clearsupcases,3),'r+');
title('clear features Levy'); view([190, -35, 110]);

% concatenate features matrix
feat_Levy = [pdf_slopes fp_binary disp_slopes];%  ];
% compute cross-correlation for exploration
CC = corr(feat_Levy,feat_Levy);
% compute normalization factor - not necessary with pca() options below
% mean
mean_featLevy = repmat(mean(feat_Levy),size(feat_Levy,1),1);
% variance
var_featLevy = repmat(var(feat_Levy),size(feat_Levy,1),1);
norm_featLevy = (feat_Levy-mean_featLevy)./var_featLevy;
% quick PCA analysis
%[wcoeff,score,latent,tsquared,explained] = pca(features_Levy,'VariableWeights',1./nfact_featLevy);
%[wcoeff,score,latent,tsquared,explained] = pca(feat_Levy,'Centered','off');
[wcoeff,score,latent,tsquared,explained] = pca(feat_Levy,'Centered','on','VariableWeights','variance');

figure; subplot(1,2,1);
plot3(score(:,1),score(:,2),score(:,3),'k+'); 
hold on; plot3(score(clearsup,1),score(clearsup,2),score(clearsup,3),'r+');
plot3(score(clearsub,1),score(clearsub,2),score(clearsub,3),'+');
xlabel('pca1'); ylabel('pca2'); zlabel('pca3');
title('features Levy'); view([190, -35, 110]);

subplot(1,2,2); plot(wcoeff(:,1:3));
ylabel('wcoeff');
xlabel('time (features)');

figure; plot(score(:,1),score(:,2),'k+'); 
hold on; plot(score(clearsup,1),score(clearsup,2),'r+');
plot(score(clearsub,1),score(clearsub,2),'+');
