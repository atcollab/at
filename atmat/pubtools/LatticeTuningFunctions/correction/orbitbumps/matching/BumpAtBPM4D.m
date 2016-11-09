function [rbump,hs,vs]=BumpAtBPM4D(ring0,inCOD,bumph,bumpv,indBPMbump,indHCor,indVCor,doplot)
% function roff=BumpAtBPM(...
%     ring0,...     1) AT lattice structure
%     inCOD,...     2) initial 6x1 coordinate guess (unused)
%     bumph,...     3) hor. bump value at indBPMbump
%     bumpv,...     4) ver. bump value at indBPMbump
%     indBPMbump,...5) bump position
%     indHCor,....  6) 3x1 correctors to generate bump. last is used for COD=0
%     indVCor,....  7) 3x1 correctors to generate bump. last is used for COD=0
%     doplot)      %8) output figure
%
% ex:
%      % order of correctors does not metter as far as the bpm is within
%      the three correctors. last corrector index is used to match the
%      postion and angle back to zero
%      roff=BumpAtBPM(ring0,0.0,1e-3,50,[4 78 90],[89 34 1]);
%
%      % to match bump at first bpm, use last corrector,
%      roff=BumpAtBPM(ring0,1e-3,1e-7,1,indHCor([end,1,2]),indVCor([end,1,2]));
% 
%see also: atmatch
if nargin<8
    doplot=false;
end

h1=atVariableBuilder(ring0,indHCor(1),{'PolynomB',{1,1}});
h2=atVariableBuilder(ring0,indHCor(2),{'PolynomB',{1,1}});
h3=atVariableBuilder(ring0,indHCor(3),{'PolynomB',{1,1}});
v1=atVariableBuilder(ring0,indVCor(1),{'PolynomA',{1,1}});
v2=atVariableBuilder(ring0,indVCor(2),{'PolynomA',{1,1}});
v3=atVariableBuilder(ring0,indVCor(3),{'PolynomA',{1,1}});
VariabH=[h1 h2 h3];
VariabV=[v1 v2 v3];

% 4D orbit
LinConstr1h=atlinconstraint(...
    indBPMbump,...
    {{'ClosedOrbit',{1}}},...
    bumph,...
    bumph,...
    1e-6);

LinConstr2h=atlinconstraint(...
    indHCor(end)+1,...
    {{'ClosedOrbit',{1}},{'ClosedOrbit',{2}}},...
    [0,0],...
    [0,0],...
    [1e-6 1e-6]);

LinConstr1v=atlinconstraint(...
    indBPMbump,...
    {{'ClosedOrbit',{3}}},...
    bumpv,...
    bumpv,...
    1e-6);

LinConstr2v=atlinconstraint(...
    indVCor(end)+1,...
    {{'ClosedOrbit',{3}},{'ClosedOrbit',{4}}},...
    [0,0],...
    [0,0],...
    [1e-6 1e-6]);
ConstrH=[LinConstr1h,LinConstr2h];
ConstrV=[LinConstr1v,LinConstr2v];


rbump=ring0;
rbump=atmatch(rbump,VariabH,ConstrH,10^-10,50,3,@lsqnonlin);%,'fminsearch');%
rbump=atmatch(rbump,VariabV,ConstrV,10^-10,50,0,@lsqnonlin);%,'fminsearch');%

if doplot
ib=findcells(ring0,'Class','Monitor');

    figure;
    o0=findorbit4Err(ring0,0,ib);
    o=findorbit4Err(rbump,0,ib);
    s=findspos(ring0,ib);
    plot(s,(o-o0)');
    legend('x','xp','y','yp');
    xlabel('s [m]');
    ylabel('4D coordinates')
    
end

end
