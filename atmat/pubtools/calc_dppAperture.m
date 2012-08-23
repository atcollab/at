function [deltap,deltam, indextab] = calc_dppAperture(THERING, varargin)
%calculate the momentum aperture at each location of the ring due to
%transverse dynamic. 
%Tracking including longitudinal motion and radiation damping
%[deltap,deltam, indextab] = calc_dppAperture.m(THERING)
%[deltap,deltam, indextab] = calc_dppAperture.m(THERING, PAmarker, PAlim)
% For SPEAR PAmarkder = 'SEPTUM'; PAlim=[-0.025 0.025 -0.005 0.005]
%[deltap,deltam, indextab] = calc_dppAperture.m(THERING, PAmarker, PAlim, Vrf)
% Vrf in MV
%output:
%   deltap, deltam, lists of postive and negative momentum aperture at
%   elements indexed by indextab
%

if nargin>=3
   PAmarkder = varargin{1};
   PAlim = varargin{2};

   ati = atindex(THERING);
   sepindx = eval(['ati.' PAmarkder]); %ati.SEPTUM;
   THERING{sepindx}.Limits = PAlim; %[-0.025 0.025 -0.005 0.005];   %apply physical aperture
   THERING{sepindx}.PassMethod = 'AperturePass';
end

ati = atindex(THERING);
Vrf = 3200000;
if nargin>=4
   Vrf = varargin{3}*1.0e6; %Volt
end
%THERING = setcellstruct(THERING,'Voltage', ati.RF, Vrf); 
%THERING = setcellstruct(THERING,'PassMethod', ati.RF, 'ThinCavityPass'); 
%turn cavity on
%cavityon(3e9);

%turn radiation on
%bi = findcells(THERING,'PassMethod','BendLinearPass');
%bi = sort([ati.BND, ati.B34]);
bi = findcells(THERING,'BendingAngle');
k = getcellstruct(THERING,'K',bi);
THERING = setcellstruct(THERING,'Energy', bi, 3e9); 
THERING = setcellstruct(THERING,'PolynomA',bi,0,1,2);
THERING = setcellstruct(THERING,'PolynomB',bi,k,1,2);
THERING = setcellstruct(THERING,'NumIntSteps',bi,10);
THERING = setcellstruct(THERING,'MaxOrder',bi,1);
THERING = setcellstruct(THERING,'PassMethod',bi,'BndMPoleSymplectic4RadPass');

%allspos = findspos(THERING, 1:length(THERING));
[curH,allspos]= getcurH(THERING);
curH(length(allspos)+1)=curH(length(allspos)); %add an extra element for the last value.
cnt = 0;
RING0 = THERING;
epsH=1e-7;
tic
for ii=1:length(allspos)
   if (RING0{ii}.Length>0 && abs(curH(ii+1)-curH(ii))>epsH)
      
       cnt = cnt + 1;
      
      indextab(cnt) = ii;
      THERING = [RING0(ii:end), RING0(1:ii-1)];
      [dplus,dminus] = track4dppAp(THERING, 500)
      deltap(cnt) = dplus;
      deltam(cnt) = dminus;
   end
end
toc

[deltap,deltam] = atmomap(ring, Nturn);
