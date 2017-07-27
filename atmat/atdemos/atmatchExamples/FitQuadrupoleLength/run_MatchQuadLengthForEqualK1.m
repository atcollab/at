% run the matching on S10 ESRF upgrade arc lattice

load('dba.mat','RING');
addpath(fullfile(pwd,'..'))
%arc=[RING(1:18);RING(128:end)];
arc0=RING;

DK1_b=compK1(arc0,indQFM,indQF);
LQ6_b=getcellstruct(arc0,'Length',indQF)*length(indQF);
LQ8_b=getcellstruct(arc0,'Length',indQFM)*length(indQFM);

% slowly reduce gradient gap, changing length of 2 quadrupoles.
% during the optimization the cell optics are kept constant retuing all
% quadrupoles. 
arc1=arc0;
warning off;
for deltaKQ=[(DK1_b-0.1):-0.2:0.0 0.0]
arc1=MatchQuadLengthForEqualK1(arc1,deltaKQ);
end
warning on;

indQF=findcells(arc1,'FamName','QF') ;
indQFM=findcells(arc1,'FamName','QFM') ;


DK1_a=compK1(arc1,indQF,indQFM);
LQ6_a=getcellstruct(arc1,'Length',indQF)*length(indQF);
LQ8_a=getcellstruct(arc1,'Length',indQFM)*length(indQFM);


disp('      before     after')
disp([['DK1  ';'LQF  ';'LQFM ';'sum L'],...
    num2str([[DK1_b;LQ6_b(1);LQ8_b(1);LQ6_b(1)+LQ8_b(1)],...
    [DK1_a;LQ6_a(1);LQ8_a(1);LQ6_a(1)+LQ8_a(1)]],'   %1.5f')])

