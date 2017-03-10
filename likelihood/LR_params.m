insampfrac = 1;
crossval = 'testset';
numchiefdomSLE = 153;

prunetime = 0;

filtertime = 3000;

gravity_type = 'input'; % 'input', 'classical' or 'grav2param' or 'grav3param'
levyfit = 'fixalpha'; % 'fixalpha' 'clauset' 'modelfun'

alevy = 2.1; %linspace(1,3,21);
rho = 2.1;
tau1 = 1;
tau2 = 0.1;

resting = 'option3'; 
chiefpop = 'median';

xcutvalue = 1;
xminvalue = 50;
xmaxvalue = 500;

rad_chiefdom = 10;

timertype = 'window'; % 'cumulative' or 'window'
cumultpts = 1;
twinsize = 60;
dayslide = 60;
datebounddefault = 0;

if datebounddefault==1
    dateboundlow =  735676;%min(ebolaTime);
    dateboundhigh = 736220;%max(ebolaTime);
else
    dateboundlow = 735676+300;
    dateboundhigh = 735676+360;
end

xminkm = 1;
% this is roughly the average radius of a chiefdom if they are round and
% equally-sized, which they are not
