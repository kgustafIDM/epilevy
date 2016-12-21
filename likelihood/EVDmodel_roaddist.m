%%
addpath C:\Users\kgustafson\EVD\Rambaut\space-time-master\space-time-master
addpath U:\MATLAB\plbinnedcode
addpath(genpath('C:\Users\kgustafson\ihv012supp_data_FamulareHu\'))
addpath U:\EpiLevy\code_base\epilevy\
addpath U:\EpiLevy
addpath C:\Users\kgustafson\EVD\kyle
addpath(genpath('C:\Users\kgustafson\Tajikistan'))
load('C:\Users\kgustafson\EVD\Rambaut\space-time-master\space-time-master\8Sept2016.mat', 'ebolaAdminL2')
load('C:\Users\kgustafson\EVD\Rambaut\space-time-master\space-time-master\8Sept2016.mat', 'ebolaAdminL3')
%%

plotLRon = 1;
needtoload = 0;

phylotype = 'tree'
dataset = 'Rambaut'

alevy = 2.1; %linspace(1,3,21);
rho = 0.0;
tau1 = 0;
tau2 = 0;
gravity_type = 'ínput'; % ínput', 'classical' or 'grav2param' or 'grav3param'

levyfit = 'fixalpha'

chiefpop = 'median'
timertype = 'window' 

cumultpts = 1;
twinsize = 30;
dayslide = 30;
datebounddefault = 1;

if datebounddefault==1
    dateboundlow = 735676;%min(ebolaTime);
    dateboundhigh = 736220;%max(ebolaTime);
else
    dateboundlow = 735676+210-30;
    dateboundhigh = 735676+210+30;
end

numchiefdomSLE = 153;

xminkm = 1;
mindd = 50; maxdd = 500;
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
%%
gg0 = 0.1:0.5:2.1;
gp1 = 0.15:0.1:1.35;
gp2 = 0.15:0.1:1.35;

plotLRon = 0;

levyfit = 'fixalpha';

chiefpop = 'median'
timertype = 'window' 

twinsize = 60;
dayslide = 60;

gravity_type = 'classical'; % 'classical' or 'grav2param' or 'grav3param'

LR_linkage = zeros(1,3+9);
LR_linkagenorm = zeros(1,3+9);
pscore_link = zeros(1,3+9);

for sss=1:size(gg0,2)
    for ttt=1:size(gp1,2)
        for uuu=1:size(gp2,2)

            rho = gg0(sss);
            tau1 = gp1(ttt);
            tau2 = gp2(uuu);
            
            alevy = 2.1;
            
            LR_linkage_plot
            
            LR_linkagenorm = [LR_linkagenorm; [rho tau1 tau2 normRset]];
            LR_linkage = [LR_linkage; [rho tau1 tau2 RL1L2set]];
            pscore_link = [pscore_link; [rho tau1 tau2 Pset]];
        end
    end
end
%%
gg0 = 2.0;
gp1 = 1.0;
gp2 = 1.0;
% gg0 = 2;
% gp1 = 1;
% gp2 = 1;
% gg0 = 2.1;
% gp1 = 0.;
% gp2 = 0.;

sss = 1;
rho = gg0(sss);
tau1 = gp1(sss);
tau2 = gp2(sss);

alevy = 2.1;

plotLRon = 1;

levyfit = 'fixalpha'; % 'fixalpha', 'clauset' 'modelfun'

chiefpop = 'median';
timertype = 'window';

gravity_type = 'input'; % 'input', 'classical' or 'grav2param' or 'grav3param'

twinsize = 60;
dayslide = 60;

dateboundlow = 735676+210-30;
dateboundhigh = 735676+210+30;
% dateboundlow = 735676;%min(ebolaTime);
% dateboundhigh = 736220;%max(ebolaTime);

LR_linkage_plot
