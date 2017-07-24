function undulator=atundulator(LUnd,nperiod,varargin)
% define undulator model 
% 
% input:
%   Lund= undulator length
%   nperiod = number of periods
%   'BendAngle', value : half pole bending angle in rad
%   'B0andEnergy', value (2x1): [half pole B0 field in T, Energy in eV]
%                               converts to bending angle in rad.
%   'magnetmodel', value : 'multipoles' (default) or 'rectangularbend'
%   'PoleGap', value : drift space between gaps defualt (0.0) 
%
% if neither BendAngle nor B0andEnergy are provided then 'BendAngle' is 0.0 
% 
% output: 
%    cellarray of elements describing un undulator of length LUnd divided 
%    in nperiod periods, each described as follows: 
%    [negpole,drift,pospole,pospole,drift,negpole]
%    if 'PoleGap' is 0.0 (default), then
%    [negpole,pospole,pospole,negpole]
%
% example: 
% 1)  und=atundulator(1.6,61,'B0andEnergy',[0.4 6.04e9])
% 2)  und=atundulator(1.6,61,'BendAngle',-0.007984472464733)
% 3)  und=atundulator(1.6,61,'B0andEnergy',[0.4 6.04e9],'magnetmodel','rectangularbend')
% 4)  und=atundulator(1.6,61,'B0andEnergy',[0.4 6.04e9],'PoleGap',0.001);
%
%see also: 

defaultAngPole=NaN;
defaultB0andEnergy=[NaN NaN];
defaultPoleGap=0;

expectedmagmodels={'multipoles','rectangularbend'};

p= inputParser;
addRequired(p,'LUnd',@isnumeric);
addRequired(p,'nperiod',@isnumeric);
addParameter(p,'magnetmodel',expectedmagmodels{1},@(x)any(validatestring(x,expectedmagmodels)));
addParameter(p,'BendAngle',defaultAngPole,@isnumeric);
addParameter(p,'B0andEnergy',defaultB0andEnergy,@isnumeric);
addParameter(p,'PoleGap',defaultPoleGap,@isnumeric);

parse(p,LUnd,nperiod,varargin{:});
magnetmodel=p.Results.magnetmodel;
BendAngle=p.Results.BendAngle;
B0andEnergy=p.Results.B0andEnergy;
PoleGap=p.Results.PoleGap;

% length of one period
periodL=LUnd/nperiod;

DistPole=PoleGap;
LPole=(periodL-2*DistPole)/4;

if ~isnan(B0andEnergy)
    B=B0andEnergy(1);%0.6;
    Brho=B0andEnergy(2)/299792458;
    AngPole=(B*LPole)/Brho;
elseif ~isnan(BendAngle)
    AngPole=BendAngle;
else
    AngPole=0.0;
end

switch magnetmodel
    
    case 'multipoles'
       
        undperiod=makeundperiod(...
            atmultipole('NegPole',LPole,0,-AngPole/LPole),...
            atmultipole('PosPole',LPole,0,AngPole/LPole),...
            atdrift('PoleGap',DistPole));
        
    case 'rectangularbend'
        
        undperiod=makeundperiod(...
            atrbend('NegPole',LPole,-AngPole,0,'BndMPoleSymplectic4Pass'),...
            atrbend('PosPole',LPole,AngPole,0,'BndMPoleSymplectic4Pass'),...
            atdrift('PoleGap',DistPole));
        
end

undulator=repmat(undperiod,nperiod,1);

end


function undper=makeundperiod(halfnegpole,halfpospole,driftpole)

if driftpole.Length>0
    undper={...
        halfnegpole;...
        driftpole;...
        halfpospole;...
        halfpospole;...
        driftpole;...
        halfnegpole;...
        };
elseif driftpole.Length==0
    undper={...
        halfnegpole;...
        halfpospole;...
        halfpospole;...
        halfnegpole;...
        };
end

end
        