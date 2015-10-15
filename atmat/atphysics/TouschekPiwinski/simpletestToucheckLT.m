% simple example of use of toucheck lifetime formula:
% [Tl,contributionsTL]=TouschekPiwinskiLifeTime(r,dpp,Ib,varargin)

close all

global GLOBVAL
GLOBVAL.E0=0.51e9;
E0=0.51e9;

load('dba.mat','RING'); %small DBA ring

%% add rf cavity 
rf.('HarmNumber')=4;%*32
rf.('Length')=0;%*32
rf.('Energy')=E0;%*32
rf.('PassMethod')='IdentityPass';%*32
rf.('Class')='Cavity';%*32
cl=299792458;
larc=findspos(RING,length(RING)+1);
rf.('Frequency')=cl/larc*rf.('HarmNumber');
rf.('Voltage')=0.5e6;%rf.('Voltage');
RING=[{rf};RING];


%% define relevant positions for evaluation of momentum aperture and LT
% this positions are default selection in the TouschekPiwinskiLifeTime.m
% function
positions=findcells(RING,'Length');
L=getcellstruct(RING,'Length',positions);
positions=positions(L>0);
    
%% get energy sperad and bunch length
RING=atradon(RING);
[l,a]=atx(RING,0,1);
%a.espread;
%a.blength;

nturns=5; %% put 500! 5 is only for speed reason in the test
%% get momentum aperture 
%(this or any better function)
[dppM,dppP]=MomAperture_allRing(RING,positions,nturns);%

% current per bunch

Ib=0.2; % mA

%% evaluate Life Time

% basic mode, constant momentum aperture
TL=TouschekPiwinskiLifeTime(RING,0.03,Ib)/3600;
disp(' basic mode, constant momentum aperture: ')
disp(TL)
% one sided momentum aperture.
TL=TouschekPiwinskiLifeTime(RING,dppM,Ib)/3600;
disp('one sided momentum aperture :')      
disp(TL)      
% two sided momentum aperture: 1/Ttot=1/2(1/Tup+1/Tdown)
TL=TouschekPiwinskiLifeTime(RING,[dppM dppP],Ib)/3600;
disp('two sided momentum aperture: 1/Ttot=1/2(1/Tup+1/Tdown) :')      
disp(TL)
% input positions
TL=TouschekPiwinskiLifeTime(RING,[dppM dppP],Ib,positions)/3600;
disp('input positions :')      
disp(TL)
% input positions and emittances
TL=TouschekPiwinskiLifeTime(RING,[dppM dppP],Ib,positions,10e-9,10e-12)/3600;
disp('input positions and emittances :')      
disp(TL)
% input positions and emittances and integration method
TL=TouschekPiwinskiLifeTime(RING,[dppM dppP],Ib,positions,10e-9,10e-12,'quad')/3600;
disp('input positions and emittances and integration method :')      
disp(TL)
% input positions and emittances and integration method and bumchelenght
% and energy spread
TL=TouschekPiwinskiLifeTime(RING,[dppM dppP],Ib,positions,10e-9,10e-12,'quad',a.espread,a.blength)/3600;
disp('input positions and emittances and integration method and bunch lenght and energy spread :')      
disp(TL)

