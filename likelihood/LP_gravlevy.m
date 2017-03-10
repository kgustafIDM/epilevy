function [logLgrav,logLlevy,logLRL1L2,normLR12,p12,rprob,mprob,gravgamma,gravprobnorm,levyprobnorm] = ...
    LP_gravlevy(rho,tau1,tau2,alevy,fgcoeff,fgcoef2,dd_meters_fullsym,cases,oL3,dL3,SLEpop,nL3SLE,snL3SLE,aX3cdfmore,aX3cdfless)

LR_params;

% allow flexibility in defintion of driving distance
drive_dist = dd_meters_fullsym./1000;

if cases>0
    if strcmp(gravity_type,'grav3param')
        gravp1 = fgcoeff.Estimate(2);
        gravp2 = fgcoeff.Estimate(3);
        gravgamma = fgcoeff.Estimate(4);
    elseif strcmp(gravity_type,'grav2param')
        gravp1 = fgcoef2.Estimate(2);
        gravp2 = fgcoef2.Estimate(3);
        gravgamma = rho;
        %     gravgamma = mdlpow.Coefficients.Estimate(2);
    elseif strcmp(gravity_type,'classical')
        gravgamma = 2; gravp1 = 1; gravp2 = 1;
    elseif strcmp(gravity_type,'input')
        gravgamma = rho; gravp1 = tau1; gravp2 = tau2;
    end
else
    gravp1 = 0;
    gravp2 = 0;
    gravgamma = 0;
end

sigmaRgravlevy = zeros(size(alevy));
p12 = zeros(size(alevy));
normLR12 = zeros(size(alevy));

logLRL1L2 = zeros(size(alevy)); pctdffR = zeros(size(alevy));
clear originid sortuniqorigins destys
% indices from the linkage data to reference the driving distance matrix
if strcmp(chiefpop,'max')
    originid = oL3.max;
    destinationid = dL3.max;
elseif strcmp(chiefpop,'min')
    originid = oL3.min;
    destinationid = dL3.min;
elseif strcmp(chiefpop,'mean')
    originid = oL3.mean;
    destinationid = dL3.mean;
elseif strcmp(chiefpop,'median')
    originid = oL3.median;
    destinationid = dL3.median;
else
    error('chiefpop not specified properly');
end
% only compute for each origin once
[sortuniqorigins,isort,iorig] = unique(originid);

rprob = zeros(numchiefdomSLE,1);
mprob = zeros(numchiefdomSLE,1);
rests = zeros(size(sortuniqorigins));
moves = zeros(size(sortuniqorigins));
for ii = 1:size(sortuniqorigins,1)
    % recover the chiefdom id for the origin
    oid = sortuniqorigins(ii);
    destys = destinationid(find(oid==originid));
    for jj = 1:size(destys)
        if oid == destys(jj)
            rests(ii) = rests(ii) + 1;
        else
            moves(ii) = moves(ii) + 1;
        end
    end
end

rprob(sortuniqorigins) = rests./(rests+moves);
mprob(sortuniqorigins) = moves./(rests+moves);
%
% restprob = restprob./100000;
% moveprob = moveprob./100000;

pop_counts = zeros(size(SLEpop,1),1);
for ii = 1:size(SLEpop,1)
    pop_counts(ii) = SLEpop(find(strcmp(nL3SLE(ii),snL3SLE)));
end

for kk = 1:size(alevy,2)
    
    if strcmp(levyfit,'fixalpha')
        alphalevy = alevy(kk);
        alphalevyless = alphalevy;
        alphalevymore = alphalevy;
        %         alphalevy = gravgamma;
    elseif strcmp(levyfit,'clauset')
        alphalevymore = aX3cdfmore;
        alphalevyless = aX3cdfless;
    elseif strcmp(levyfit,'modelfun')
        alphalevy = mdlpow.Coefficients.Estimate(2);
    end
    
    gravprob = zeros(size(drive_dist));
    levyprob = zeros(size(drive_dist));
    
    gravprobnorm = zeros(size(drive_dist));
    levyprobnorm = zeros(size(drive_dist));
    
    normgrav = zeros(size(drive_dist,1),1);
    normlevy = zeros(size(drive_dist,1),1);
    
    for ii = 1:size(drive_dist,1)
        for jj = 1:size(drive_dist,1)
            if ii == jj
                if strcmp(resting,'none')
                    gravprob(ii,jj) = 0;
                    levyprob(ii,jj) = 0;
                elseif strcmp(resting,'option1')
                    gravprob(ii,jj) = pop_counts(ii)^2./rad_chiefdom^2;
                    levyprob(ii,jj) = 1/rad_chiefdom^alphalevy;
                elseif strcmp(resting,'option2')
                    %                     gravprob(ii,jj) = 1;
                    %                     gravprob(ii,jj) = pop_counts(ii)^(2*gravp1)./rad_chiefdom^gravgamma;
                    %                     levyprob(ii,jj) = 1/rad_chiefdom^alphalevy;
                    gravprob(ii,jj) = pop_counts(ii);
                    levyprob(ii,jj) = 1/rad_chiefdom^alphalevy;
                    %                     levyprob(ii,jj) = mean(restprob);
                    %                     levyprob(ii,jj) = pop_counts(ii)./rad_chiefdom;
                end
            else
                gravprob(ii,jj) = pop_counts(ii)^gravp1*pop_counts(jj)^gravp2/(drive_dist(ii,jj)/xminkm)^gravgamma;
                % experimenting with assortative Levy flight model
                if(drive_dist(ii,jj)/xminkm < xcutvalue)
                    levyprob(ii,jj) = 1/(drive_dist(ii,jj)/xminkm)^alphalevyless;
                elseif(drive_dist(ii,jj)/xminkm > xcutvalue)
                    levyprob(ii,jj) = 1/(drive_dist(ii,jj)/xminkm)^alphalevymore;
                end
            end
        end
        % this is the chosen option for publication
        if strcmp(resting,'option3')
            gravprob(ii,ii) = sum(gravprob(ii,:)); % forces resting probability to 1/2 the total probability
            levyprob(ii,ii) = sum(levyprob(ii,:)); % forces resting probability to 1/2 the total probability
        end
        normgrav(ii) = sum(gravprob(ii,:));
        normlevy(ii) = sum(levyprob(ii,:));
        gravprobnorm(ii,:) = gravprob(ii,:)./normgrav(ii);
        levyprobnorm(ii,:) = levyprob(ii,:)./normlevy(ii);
        % make the (ii,ii) resting probablity for Levy equal to gravity,
        % and rescale the rest of the Levy probabilities to conserve total
        %         levyprobnorm(ii,ii) = gravprobnorm(ii,ii);
        %         sumelse = sum(levyprobnorm(ii,[1:ii-1 ii+1:end]));
        %         levyprobnorm(ii,[1:ii-1 ii+1:end]) = levyprobnorm(ii,[1:ii-1 ii+1:end])*(1-levyprobnorm(ii,ii))./sumelse;
    end
    
    %compute model likelihoods
    % start with blank vector -  then fill with every value for
    % origin-destination pairs
    
    logLgrav = [];
    logLlevy = [];
    
    opopsy = [];
    for ii = 1:size(sortuniqorigins,1)
        % recover the chiefdom id for the origin
        oid = sortuniqorigins(ii);
        opopsy(ii) = pop_counts(oid);
        destys = destinationid(find(oid==originid));
        for jj = 1:size(destys,1)
            ddij = drive_dist(oid,destys(jj));
            if oid == destys(jj)
                % this protects from dividing by zero
                ijlogLlevy = log(gravprobnorm(oid,oid));
                ijlogLgrav = log(levyprobnorm(oid,oid));
            else
                % evaluate the normalized gravity probability at the observed
                % origin and destination from the POTN
                normedgrav = gravprobnorm(oid,destys(jj));
                % take the logarithm
                ijlogLgrav =  log(normedgrav);
                % evaluate the normalized levy probability at the observed
                % origin and destination from the POTN
                clear nextterm1 nextterm2
                if(ddij > xcutvalue) % implement assortative model
                    nextterm1 = -log(normlevy(oid));
                    nextterm2 = -alphalevymore*log(ddij/xminkm);
                    ijlogLlevy = nextterm1 + nextterm2;
                elseif(ddij < xcutvalue) % set boundary values for Levy probability
                    nextterm1 = -log(normlevy(oid));
                    nextterm2 = -alphalevyless*log(ddij/xminkm);
                    ijlogLlevy = nextterm1 + nextterm2;
                end
            end
            logLgrav = [logLgrav ijlogLgrav];
            logLlevy = [logLlevy ijlogLlevy];
        end
    end
    
    l1 = logLgrav;
    l2 = logLlevy;
    
    % Calculate likelihood ratio (LR)
    % since this is grav - Levy, a positive value prefers gravity
    logLRL1L2(kk) = sum(l1)-sum(l2);
    pctdffR(kk) = 100.*(sum(l1)-sum(l2))/sum(l1);
    
    nsize = size(l1,2);
    
    % Calculate mean of log likelihoods
    l1_bar = (1/nsize)*sum(l1);
    l2_bar = (1/nsize)*sum(l2);
    
    temp = ((l1-l2)-(l1_bar - l2_bar)).^2;
    sigmaRgravlevy(kk) = sqrt((1/nsize)*sum(temp(~isnan(temp))));
    % Vuoung test
    p12(kk) = erfc(abs(logLRL1L2(kk))/(sqrt(2*nsize)*sigmaRgravlevy(kk)));
    if p12(kk)==0
        % minimum possible value of the erfc computable to five sigfigs
        p12(kk) =  4.9407e-324;
    end;
    if abs(logLRL1L2) < 1E-10
        normLR12(kk)=0;
    else
        normLR12(kk) = logLRL1L2(kk)/(sqrt(nsize)*sigmaRgravlevy(kk));
    end
    
end

% figure; subplot(1,2,1);
% plot(alevy,RL1L2); xlabel('Levy \alpha'); ylabel('Ratio'); grid on
% subplot(1,2,2); semilogy(alevy,pgravlevy); xlabel('Levy \alpha'); ylabel('p-value');
if size(alevy,2)>1
    figure; plotyy(alevy,logLRL1L2,alevy,-log10(p12)); grid on;
    xlabel('\alpha for Levy flight');
else
    logLRL1L2;
    -log10(p12);
end

end