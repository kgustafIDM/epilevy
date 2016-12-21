%%

% allow flexibility in defintion of driving distance
drive_dist = drive_dist_meters_fullsym./1000;

if casesat(k)>0 && strcmp(gravity_type,'grav3param')
    gravp1 = fitgravcoeff{1, k}.Estimate(2);
    gravp2 = fitgravcoeff{1, k}.Estimate(3);
    gravgamma = fitgravcoeff{1, k}.Estimate(4);
elseif casesat(k)>0 && strcmp(gravity_type,'grav2param')
    gravp1 = fitgravcoef2{1, k}.Estimate(2);
    gravp2 = fitgravcoef2{1, k}.Estimate(3);
    gravgamma = 2;
%     gravgamma = mdlpow.Coefficients.Estimate(2);
elseif strcmp(gravity_type,'classical')
    gravgamma = 2; gravp1 = 1; gravp2 = 1;
elseif strcmp(gravity_type,'input')
    gravgamma = rho; gravp1 = tau1; gravp2 = tau2;
end

sigmaRgravlevy = zeros(size(alevy));
pgravlevy = zeros(size(alevy));
normRgravlevy = zeros(size(alevy));

RL1L2 = zeros(size(alevy)); pctdffR = zeros(size(alevy));
clear originid sortuniqorigins destys
% indices from the linkage data to reference the driving distance matrix
if strcmp(chiefpop,'max')
    originid = originL3_max;
    destinationid = destinationL3_max;
elseif strcmp(chiefpop,'min')
    originid = originL3_min;
    destinationid = destinationL3_min;
elseif strcmp(chiefpop,'mean')
    originid = originL3_mean;
    destinationid = destinationL3_mean;
elseif strcmp(chiefpop,'median')
    originid = originL3_median;
    destinationid = destinationL3_median;
else
    error('chiefpop not specified properly');
end
% only compute for each origin once
[sortuniqorigins,isort,iorig] = unique(originid);

restprob = zeros(numchiefdomSLE,1);
moveprob = zeros(numchiefdomSLE,1);
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

restprob(sortuniqorigins) = rests./(rests+moves);
moveprob(sortuniqorigins) = moves./(rests+moves);
% 
% restprob = restprob./100000;
% moveprob = moveprob./100000;

pop_counts = zeros(size(SLPpop,1),1);
for ii = 1:size(SLPpop,1)
    pop_counts(ii) = SLPpop(find(strcmp(nameL3SLE(ii),sortednameL3SLE)));
end

for kk = 1:size(alevy,2)
    
    if strcmp(levyfit,'fixalpha')
        alphalevy = alevy(kk);
%         alphalevy = gravgamma;
    elseif strcmp(levyfit,'clauset')
        alphalevy = aX3cdf;
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
                elseif strcmp(resting,'option1');
%                     gravprob(ii,jj) = pop_counts(ii)^(2*gravp1)./rad_chiefdom^gravgamma;
%                     gravprob(ii,jj) = pop_counts(ii)^2./rad_chiefdom^2;
%                     gravprob(ii,jj) = 1;
%                     gravprob(ii,jj) = pop_counts(ii)./rad_chiefdom;

                    gravprob(ii,jj) = 0;
                    levyprob(ii,jj) = 0;

%                     levyprob(ii,jj) = restprob(ii);
%                     levyprob(ii,jj) = mean(restprob);
%                     levyprob(ii,jj) = 1/rad_chiefdom^alphalevy;
%                     levyprob(ii,jj) = pop_counts(ii)./rad_chiefdom;
                end
            else
                gravprob(ii,jj) = pop_counts(ii)^gravp1*pop_counts(jj)^gravp2/(drive_dist(ii,jj)/xminkm)^gravgamma;
                levyprob(ii,jj) = 1/(drive_dist(ii,jj)/xminkm)^alphalevy;
            end
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
    normg = [];
    pgrav = [];
    normlev = [];
    plev = [];
    
    opopsy = [];
    for ii = 1:size(sortuniqorigins,1)
        % recover the chiefdom id for the origin
        oid = sortuniqorigins(ii);
        opopsy(ii) = pop_counts(oid);
        destys = destinationid(find(oid==originid));
        for jj = 1:size(destys,1)
            ddij = drive_dist(oid,destys(jj));
            if (ddij > 0 && ddij < maxdd)
                if oid == destys(jj)
                    % this protects from dividing by zero
                    ijlogLlevy = log(gravprobnorm(ii,ii));
                    ijlogLgrav = log(levyprobnorm(ii,ii));
                else
                    
                    normg1 = -log(normgrav(oid));
                    pgrav1 = gravp1*log(pop_counts(oid)) ...
                             + gravp2*log(pop_counts(destys(jj)))...
                             - gravgamma*log(ddij/xminkm);
                   
                    ijlogLgrav =  normg1 + pgrav1;
                    
                    nextterm1 = -log(normlevy(oid));
                    nextterm2 = -alphalevy*log(ddij/xminkm);
                    
                    ijlogLlevy = nextterm1 + nextterm2;
                    
                end
                normg = [normg normg1];
                pgrav = [pgrav pgrav1];
                normlev = [normlev nextterm1];
                plev = [plev nextterm2];
                
                logLgrav = [logLgrav ijlogLgrav];
                logLlevy = [logLlevy ijlogLlevy];
            end
        end
    end
        
    l1 = logLgrav;
    l2 = logLlevy;
    
    % Calculate likelihood ratio (LR)
    % since this is grav - Levy, a positive value prefers gravity
    RL1L2(kk) = sum(l1)-sum(l2);
    pctdffR(kk) = 100.*(sum(l1)-sum(l2))/sum(l1);
    
    nsize = size(l1,2);
    
    % Calculate mean of log likelihoods
    l1_bar = (1/nsize)*sum(l1);
    l2_bar = (1/nsize)*sum(l2);     
    
    temp = ((l1-l2)-(l1_bar - l2_bar)).^2;
    sigmaRgravlevy(kk) = sqrt((1/nsize)*sum(temp(~isnan(temp))));
    % Vuoung test
    pgravlevy(kk) = erfc(abs(RL1L2(kk))/(sqrt(2*nsize)*sigmaRgravlevy(kk)));
    if pgravlevy(kk)==0
        % minimum possible value of the erfc computable to five sigfigs
        pgravlevy(kk) =  4.9407e-324
    end;
    normRgravlevy(kk) = RL1L2(kk)/(sqrt(nsize)*sigmaRgravlevy(kk)); 
    
end

% figure; subplot(1,2,1);
% plot(alevy,RL1L2); xlabel('Levy \alpha'); ylabel('Ratio'); grid on 
% subplot(1,2,2); semilogy(alevy,pgravlevy); xlabel('Levy \alpha'); ylabel('p-value');
if size(alevy,2)>1
    figure; plotyy(alevy,RL1L2,alevy,-log10(pgravlevy)); grid on;
    xlabel('\alpha for Levy flight');
else
    RL1L2;
    -log10(pgravlevy);
end