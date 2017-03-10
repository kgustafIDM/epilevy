%%
addpath C:\Users\kgustafson\EVD\Rambaut\space-time-master\space-time-master
addpath U:\MATLAB\plbinnedcode
addpath(genpath('C:\Users\kgustafson\ihv012supp_data_FamulareHu\'))
addpath U:\EpiLevy\code_base\epilevy\likelihood
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

cleanEboladata;

gg0 = 0; gp1 = 0; gp2 = 0; alevy=0;

LR_linkage_plot

%%
plotLRon = 0;

normRtest_multi = [];
normRval_multi = [];
Ptest_multi = [];
Pval_multi = [];

for kk = 1:1000
    
    normRtest_multi = [normRtest_multi; normRtest];
    Ptest_multi = [Ptest_multi; Ptest];
    
    normRval_multi = [normRval_multi; normRval];
    Pval_multi = [Pval_multi; Pval];
    
    LR_linkage_plot;
    
end

%%
gg0 = 1:0.1:6.0; %rho
gp1 = 1; %tau1
gp2 = 0.1:0.1:3; %tau2
alpha = 0.1:1:4.1;

plotLRon = 0;
phylotype = 'tree'
dataset = 'Rambaut'

cleanEboladata;

% LR_230_rhoalpha = zeros(3,size(gg0,2)*size(alpha,2));
LR_330_rhotau2alpha = zeros(5,size(gg0,2)*size(gp2,2));
% LR_linkage = zeros(1,3+9);
% LR_linkagenorm = zeros(1,3+9);
% pscore_link = zeros(1,3+9);
counter = 0;
for rrr=1:size(alpha,2)
    for sss=1:size(gg0,2)
        %for ttt=1:size(gp1,2)
        for uuu=1:size(gp2,2)
            counter = counter + 1;
            
            rho = gg0(sss);
            %             tau1 = gp1(ttt);
            %             tau2 = gp2(uuu);
            tau1 = gp1;
            tau2 = gp2(uuu);
            
            alevy = alpha(rrr)+1;
            
            LR_linkage_plot;
            
            LR_330_rhotau2alpha(1,counter) = rho;
            LR_330_rhotau2alpha(2,counter) = alevy-1;
            LR_330_rhotau2alpha(3,counter) = tau1;
            LR_330_rhotau2alpha(4,counter) = tau2;
            LR_330_rhotau2alpha(5,counter) = normRtest;
            %             LR_230_rhoalpha(1,counter) = rho;
            %             LR_230_rhoalpha(2,counter) = alevy - 1;
            %             LR_230_rhoalpha(3,counter) = normRtest;
            %             LR_linkagenorm = [LR_linkagenorm; [rho tau1 tau2 normRset]];
            %             LR_linkagenorm = [LR_linkagenorm; [rho tau1 tau2 normRset]];
            %             LR_linkage = [LR_linkage; [rho tau1 tau2 RL1L2set]];
            %             pscore_link = [pscore_link; [rho tau1 tau2 Pset]];
        end
        %end
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
