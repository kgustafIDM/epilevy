
% manual correction for Google Maps hiccup
latlong(166,2) = 8.135583;
latlong(166,1) = -12.715;
latlong(176,2) = 8.135583;
latlong(176,1) = -12.715;
latlong(181,2) = 8.135583;
latlong(181,1) = -12.715;

[lats,iulats,idlats] = unique(latlong(:,2));
[lons,iulons,idlons] = unique(latlong(:,1));

ulonsOrdered = latlong(iulats,1);
ulatsOrdered = latlong(iulats,2);

drive_dur_seconds = zeros(size(lats,1),size(lats,1));
drive_dist_meters = zeros(size(lats,1),size(lats,1));

for kk = 1:size(lats,1)-1;
    for jj = kk+1:size(lats,1);
        orig_coord = [num2str(ulatsOrdered(kk)) ',' num2str(ulonsOrdered(kk))];
        dest_coord = [num2str(ulatsOrdered(jj)) ',' num2str(ulonsOrdered(jj))];
        mode='driving';
        api_key = 'AIzaSyCzC5cA6s7aI8AmfRGtwfLGJVs4s_vwhkE';
        
        url = ['https://maps.googleapis.com/maps/api/distancematrix/json?origins=',orig_coord,'&destinations=',dest_coord,'&mode=',mode,'&language=en-EN&sensor=false&key=',api_key''];
        
        mapinfo = webread(url);
        
        drive_dur_seconds(kk,jj) = mapinfo.rows.elements.duration.value;
        drive_dist_meters(kk,jj) = mapinfo.rows.elements.distance.value;
    end
end
%%

alllat = ebolaAdminL3.numberAttribute(:,2);
alllon = ebolaAdminL3.numberAttribute(:,1);

% correction for Bempe
a1 = find(alllon==-12.706002808285380);
a2 = find(alllat==8.121984516370414);
alllat(a2) = 8.135583;
alllon(a1) = -12.715;

% correction for Dembelia Sinkunia
a3 = find(alllon==-11.530751501260138);
alllon(a3) = -11.4508;

% correction for Neini, placed at Fankoya epicenter
alllat(70) = 9.166667;
alllon(70) = -10.966667;
% correction for Neya, placed at Kurakoro
alllat(71) = 8.9683941;
alllon(71) = -11.2034561;
% correction for Sengbe, placed at Yogomaia (via Mapcarta)
alllat(72) = 9.5833;
alllon(72) = -11.5667;
% correction for Wara Wara Bafodia, placed nearest Google Map road lock
alllat(74) = 9.6397479;
alllon(74) = -11.7069335;
% correction for Selenga, placed nearest Kpatema
alllat(109) = 8.1197;
alllon(109) = -11.7325;
% correction for Dema, placed at Bendu Wharf, no boat directions 
alllat(115) = 7.8627;
alllon(115) = -12.8966;
% correction for Jong, placed nearby, other side of river 
alllat(117) = 7.6049181;
alllon(117) = -12.2192947;
% correction for Kwamebai Krim, placed nearby road lock 
alllat(119) = 7.370466;
alllon(119) = -11.9444446;
% correction for Mano Sakrim, placed nearby road lock 
alllat(142) = 7.1998197;
alllon(142) = -11.7580677;
% correction for Panga Krim, placed near Pujehun
alllat(144) = 7.3282722;
alllon(144) = -11.753297;

drive_dur_seconds_full = zeros(size(alllat,1),size(alllat,1));
drive_dist_meters_full = zeros(size(alllat,1),size(alllat,1));

wroptions = weboptions('Timeout',60);

for kk = 1:size(alllat,1)-1;
    for jj = kk+1:size(alllat,1);
        orig_coord = [num2str(alllat(kk)) ',' num2str(alllon(kk))];
        dest_coord = [num2str(alllat(jj)) ',' num2str(alllon(jj))];
        mode='driving';
%         api_key = 'AIzaSyDswkJadwaVzGuSkJ8Z0gzBcEekqdn7Wi4';
%        api_key = 'AIzaSyCzC5cA6s7aI8AmfRGtwfLGJVs4s_vwhkE';
%         api_key = 'AIzaSyAtXVxIbwYVClh1zowtbrDQOLcuF_yEsfY';
%          api_key = 'AIzaSyC_mdkFSXd4QLl6FAPgFuYnofkKkIXzWMI';
         api_key = 'AIzaSyBnynPNxpGj_9ltOjPkqaUBKz75KF66slI';
        url = ['https://maps.googleapis.com/maps/api/distancematrix/json?origins=',orig_coord,'&destinations=',dest_coord,'&mode=',mode,'&language=en-EN&sensor=false&key=',api_key''];
        
        mapinfo = webread(url,wroptions);
        
        drive_dur_seconds_full(kk,jj) = mapinfo.rows.elements.duration.value;
        drive_dist_meters_full(kk,jj) = mapinfo.rows.elements.distance.value;
    end
end

ddmf = drive_dist_meters_full(:);   
  
% result of [alpha, xmin, L] = plfit(ddmf(ddmf>0));
alpha = 3.05;
xmin = 167906;

%%
figure;
plot(ulonsOrdered,ulatsOrdered,'.r','MarkerSize',20) 
plot_google_map

figure;
plot(alllon,alllat,'.r','MarkerSize',20) 
plot_google_map

lat = ebolaAdminL3.numberAttribute(:,2); 
lon = ebolaAdminL3.numberAttribute(:,1);  
plot(lon,lat,'.k','MarkerSize',20) 
plot_google_map

%% attempt to do an entire distance matrix in one API call

clear orcoord
for kk = 1:size(alllat,1)-1;
    for jj = kk+1:size(alllat,1);
        orcoor(kk,:) = [num2str(alllat(kk),'%1.3f') ',' num2str(alllon(kk),'%1.3f')];
    end
end
lat = num2str(alllat,'%1.3f\n');
lon = num2str(alllon,'%1.3f\n');
commas = cell(153,1);
bars = cell(152,1);
for ii = 1:152
    bars{ii} = '|';
end
for ii = 1:153
    commas{ii} = ',';
end
alllall = [lat,char(commas),lon];
fullstr = [];
for jj = 1:152
    fullstr = [fullstr num2str(alllat(jj)) ',' num2str(alllon(jj,:)) '|'];
end

url = ['https://maps.googleapis.com/maps/api/distancematrix/json?origins=',fullstr,'&destinations=',fullstr,'&mode=',mode,'&language=en-EN&sensor=false&key=',api_key''];



