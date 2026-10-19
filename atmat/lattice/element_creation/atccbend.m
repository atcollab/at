function elem=atccbend(fname,L,angle,K,varargin)
%ATCCBEND Straight-axis combined-function bend with elegant CCBEND end maps.
% ELEM=ATCCBEND(NAME,L,ANGLE,K1,'Name',Value,...)
% L is the nominal arc length [m]; ANGLE is nonzero [rad].
% PolynomB uses AT coefficients: [0 K1 K2/2 K3/6 ...].
% FringeModel: -1=off, 0=legacy elegant, 1=Lindberg integral map.
% fringeInt1/2: [K0 K2 K4 K5 K6 K7 I0 I1], one array per end.
% NumIntSteps=4; IntegrationOrder=4 (also 2/6); EdgeOrder=3.
% OptimizeFSE=1; OptimizeDX=1; CompensateKN=0; ReferenceCorrection=1.
% FullGap and FringeInt1/2 provide the legacy finite-gap correction.
% PassMethod='CCBendRadPass' adds radiation loss; 'CCBendQuantPass' adds
% loss and Gaussian ISR. Supply Energy [eV] or use ring energy properties.
% This model does not include spin, SDDS multipole files or
% elegant-specific 3D misalignment. Use AT R1/R2/T1/T2 and apertures.
% See docs/common/ccbend.md and atintegrators/elegant.LICENSE.
if nargin<4, K=0; end
defaults={'NumIntSteps',4,'IntegrationOrder',4,'FringeModel',0,'EdgeOrder',3, ...
    'OptimizeFSE',1,'OptimizeDX',1,'CompensateKN',0,'ReferenceCorrection',1};
for i=1:2:numel(defaults)
    [defaults{i+1},varargin]=getoption(varargin,defaults{i},defaults{i+1});
end
elem=atrbend(fname,L,angle,K,'CCBendPass',defaults{:},varargin{:});
end
