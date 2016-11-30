%%
% addpath C:\Users\kgustafson\EVD\Rambaut\space-time-master\space-time-master
% addpath U:\MATLAB\plbinnedcode
% addpath(genpath('C:\Users\kgustafson\ihv012supp_data_FamulareHu\'))
% addpath U:\EpiLevy\code_base\epilevy\
% addpath(genpath('C:\Users\kgustafson\Tajikistan'))
% load('C:\Users\kgustafson\EVD\Rambaut\space-time-master\space-time-master\8Sept2016.mat', 'ebolaAdminL2')
% load('C:\Users\kgustafson\EVD\Rambaut\space-time-master\space-time-master\8Sept2016.mat', 'ebolaAdminL3')
%%

needtoload = 0;

phylotype = 'tree'
dataset = 'Rambaut'

alevy = 1.8; %linspace(1,3,21);
gravgamma0 = alevy(sss);
gravp10 = gp1(sss);
gravp20 = gp2(sss);
gravity_type = 'grav3param'; % 'classical' or 'grav2param' or 'grav3param'

levyfit = 'fixalpha'

chiefpop = 'median'
timertype = 'window' 

cumultpts = 11;
twinsize = 60;
datebounddefault = 1;

if datebounddefault==1
    dateboundlow = 735676;%min(ebolaTime);
    dateboundhigh = 736220;%max(ebolaTime);
else
    dateboundlow = 735676+210-60;
    dateboundhigh = 735676+210+60;
end

numchiefdomSLE = 153;

xminkm = 1;
mindd = 50; maxdd = 400;
% this is roughly the average radius of a chiefdom if they are round and
% equally-sized, which they are not
rad_chiefdom = 10;
resting = 'option1'; 

if needtoload
    addpath U:\EpiLevy\code_base\epilevy\
    if strcmp(dataset,'Park')
        Park_SLe_network_prelims
    elseif strcmp(dataset,'Rambaut')
        Rambaut_SLe_network_prelims
    else
        error('nodataset');
    end
end

if strcmp(phylotype,'tree')
    nodenameL2 = ebolaTree.sourceData.AdminL2;
    nodenameL3 = ebolaTree.sourceData.AdminL3;
    ebolaTime = ebolaTree.time;
    EVDlinkagesSLE = ebolaTree.network;
    EVDpvalueSLE = ebolaTree.p;
elseif strcmp(phylotype,'network')
    nodenameL2 = ebolaNetwork.sourceData.AdminL2;
    nodenameL3 = ebolaNetwork.sourceData.AdminL3;
    ebolaTime = ebolaNetwork.time;
    EVDlinkagesSLE = ebolaNetwork.network;
    EVDpvalueSLE = ebolaNetwork.p;
end

LR_linkage_plot