%%
cmap = hsv(14);
cmapg = gray(13);

[ebTsort, st2] = sort(ebolaTree.time);

origindates = ebolaTree.time(EVDlinkage(:,1));
destinationdates = ebolaTree.time(EVDlinkage(:,2));

[sortddates, sdids] = sort(destinationdates);
uniqddates = unique(sortddates);
cmap = flipud(cmap);
nodeL2sort = nodeL2id(st2);

%timeset = 1:80;
%timeset = 81:128;
%timeset = 128:size(uniqddates,1);
timeset = 1:size(uniqddates,1);
%timeset = 170:195;

mintime = 735676;
figure; hold on;
for kk = timeset
    currentid = find(destinationdates==uniqddates(kk));
    uniqddates(kk)-mintime
    for jj=1:size(currentid,1)
        XX = [origindates(currentid(jj))-mintime,destinationdates(currentid(jj))-mintime];
        YY = [15-nodeL2id(EVDlinkage(currentid(jj),1)),15-nodeL2id(EVDlinkage(currentid(jj),2))];
        plot(XX,YY,'Color',cmapg(10,:),'LineWidth',0.25);
        scatter(XX(1),YY(1),45,'MarkerEdgeColor',cmap(YY(1),:));
        scatter(XX(2),YY(2),65,'.','MarkerEdgeColor',cmap(YY(1),:));
        ax = gca; ax.YLim = [1 14]; ax.YTick = (1:14); ax.YTickLabel = [];%flipud(nameL2SLE);
    end
%     pause(0.2); drawnow;
end

% xlabel('days after 18-Mar-2014');
% title('pruned to shortest linkages');