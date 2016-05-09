% short script to find distribution of displacements
% 
% dpos_x = dpos_r(1:Np/2,:); dpos_y = dpos_r(Np/2+1:end,:);
% figure; hold on;
% for ii = 1:Np/2
%     plot(dpos_x(ii,:),dpos_y(ii,:));
% end

tgridsize = 10000;

% stand in for pos_r as a one-dimensional walk
% using pos_y here is an unsigned sequence of random variates
pos_r = pos_y;%sqrt(pos_x.^2+pos_y.^2);

firstzero = zeros(Np,1);

for ii = 1:Np
    fff = find(~pos_r(ii,:),1,'first');
    if isempty(fff)
        firstzero(ii) = size(pos_r,2);
    else
        firstzero(ii) = fff;
    end
end

% random_dir = random('uniform',0,1,size(pos_x));
% pos_x = pos_r.*sqrt(1-random_dir.^2);
% pos_y = pos_r.*random_dir;

random_angle = random('uniform',0,2*pi,size(pos_x));
pos_x = pos_r.*cos(random_angle);
pos_y = pos_r.*sin(random_angle);

% figure; 
% for ii = 1:Np
%     plot(cumsum(pos_y(ii,:)));
%     pause(0.5);
% end

trace_x = cumsum(pos_x,2);
trace_y = cumsum(pos_y,2);
%trace_r = cumsum(pos_r,2);

times1 = cumsum(time_r,2);
maxtimes = max(times1,[],2);

tgrid = linspace(0,min(maxtimes),tgridsize);

ystar = zeros(Np,length(tgrid));
xstar = zeros(Np,length(tgrid));
%rstar = zeros(Np,length(tgrid));

for ii = 1:Np
    truntimes = [0,times1(ii,1:firstzero(ii)-1)];
    truntrace = [0,trace_x(ii,1:firstzero(ii)-1)];
    xstar(ii,:) = interp1(truntimes,truntrace,tgrid);
    truntrace = [0,trace_y(ii,1:firstzero(ii)-1)];
    ystar(ii,:) = interp1(truntimes,truntrace,tgrid);
    %truntrace = [0,trace_r(ii,1:firstzero(ii)-1)];
    %rstar(ii,:) = interp1(truntimes,truntrace,tgrid);
end
clear trace_x trace_y
dpx = xstar-repmat(xstar(:,1),1,length(tgrid));
dpy = ystar-repmat(ystar(:,1),1,length(tgrid));
%dpr = rstar-repmat(rstar(:,1),1,length(tgrid));
clear xstar ystar
disp_r = sqrt(dpx.^2 + dpy.^2);
disp_x = dpx;
disp_y = dpy;
clear dpx dpy
% diffposr = diff(disp_r);
%  
% % for the entire ensemble
% stepsB = diffposr(:);
% [sup1,sup2] = hist(stepsB,50);
% 
% % on a per-walker basis
% istepsB = diffposr;

stddpos_r = std(disp_r,0,1);
stddpos_x = std(disp_x,0,1);
stddpos_y = std(disp_y,0,1);

diff_r = diff(disp_r,[],2);

walkbins = [];
num_bins = 10;
for ii=1:num_bins
    walkbins(ii,:) = 1+(ii-1)*Np/num_bins:ii*Np/num_bins;
end

%%
% figure; 
% for ii = 1:Np
%     %for jj = 1:length(tgrid)
%        %hold on;
%        plot(xstar(ii,:),ystar(ii,:),'r.');
%        set(gca,'Xlim',[-4 4],'Ylim',[-4 4]);
%        %pause(0.001);
%     %end
%     pause(0.1);
% end