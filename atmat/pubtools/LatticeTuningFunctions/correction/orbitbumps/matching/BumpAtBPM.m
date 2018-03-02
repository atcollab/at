function [rbump,hs,vs]=BumpAtBPM(ring0,inCOD,bumph,bumpv,indBPMbump,indHCor,indVCor)
% function roff=BumpAtBPM(...
%     ring0,... AT lattice structure
%     inCOD,... initial 6x1 coordinate guess 
%     bumph,... hor. bump value at indBPMbump
%     bumpv,... ver. bump value at indBPMbump
%     indBPMbump, bump position
%     indHCor,.... 1x3 correctors to generate bump. last is used for COD=0
%     indVCor.... 1x3 correctors to generate bump. last is used for COD=0
%     )
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
%see also: atmatch findorbit6Err

if size(indBPMbump)~=[1 1]
    error('indBPMbump must be size 1x1')
end
if size(indHCor)~=[1 3]
    error('indHCor must be size 1x3')
end
if size(indVCor)~=[1 3]
    error('indVCor must be size 1x3')
end

h1=atVariableBuilder(ring0,indHCor(1),{'PolynomB',{1,1}});
h2=atVariableBuilder(ring0,indHCor(2),{'PolynomB',{1,1}});
h3=atVariableBuilder(ring0,indHCor(3),{'PolynomB',{1,1}});
v1=atVariableBuilder(ring0,indVCor(1),{'PolynomA',{1,1}});
v2=atVariableBuilder(ring0,indVCor(2),{'PolynomA',{1,1}});
v3=atVariableBuilder(ring0,indVCor(3),{'PolynomA',{1,1}});
VariabH=[h1 h2 h3];
VariabV=[v1 v2 v3];

% 6D orbit
ConstrH6D=struct(...
    'Fun',@(r,~,~)get6dx(r,indBPMbump,indHCor(end)+1,inCOD),...
    'Weight',[1e-6 1e-6 1e-6],...
    'RefPoints',1,...
    'Min',[bumph 0.0 0.0],...
    'Max',[bumph 0.0 0.0]);

ConstrV6D=struct(...
    'Fun',@(r,~,~)get6dy(r,indBPMbump,indHCor(end)+1,inCOD),...
    'Weight',[1e-6 1e-6 1e-6],...
    'RefPoints',1,...
    'Min',[bumpv 0.0 0.0],...
    'Max',[bumpv 0.0 0.0]);


rbump=ring0;

try
    rbump=atmatch(rbump,VariabV,ConstrV6D,10^-16,10,3,@lsqnonlin);%,'fminsearch');%
    rbump=atmatch(rbump,VariabH,ConstrH6D,10^-16,100,3,@fminsearch);%,'fminsearch');%
catch
    rbump=atmatch(rbump,VariabH,ConstrH6D,10^-10,40,3);%,'fminsearch');%
    rbump=atmatch(rbump,VariabV,ConstrV6D,10^-10,40,3);%,'fminsearch');%
    rbump=atmatch(rbump,VariabH,ConstrH6D,10^-16,10,3,@lsqnonlin);%,'fminsearch');%
    rbump=atmatch(rbump,VariabV,ConstrV6D,10^-16,10,3,@lsqnonlin);%,'fminsearch');%
end

% plot corrector values.
hs=atgetfieldvalues(rbump,indHCor,'PolynomB',{1,1});
vs=atgetfieldvalues(rbump,indVCor,'PolynomA',{1,1});

end


function x=get6dx(r,ind1,ind2,inCOD)
o1=findorbit6Err(r,ind1,inCOD);
o2=findorbit6Err(r,ind2,inCOD);
x=[o1(1,1),o2(1,1),o2(2,1)]; % orbit at ind1, orbit and angle at ind2
end

function x=get6dy(r,ind1,ind2,inCOD)
o1=findorbit6Err(r,ind1,inCOD);
o2=findorbit6Err(r,ind2,inCOD);
x=[o1(3,1),o2(3,1),o2(4,1)]; % orbit at ind1, orbit and angle at ind2
end

