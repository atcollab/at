% run the matching on S10 ESRF upgrade arc lattice

load('dba.mat','RING');
addpath(fullfile(pwd,'..'))
%arc=[RING(1:18);RING(128:end)];
arc=RING;

% proceed in steps
arc1=MatchQuadLengthForEqualK1(arc,0.7);
arc2=MatchQuadLengthForEqualK1(arc1,0.5);
arc3=MatchQuadLengthForEqualK1(arc2,0.3);
arc4=MatchQuadLengthForEqualK1(arc3,0.15);
arc1=MatchQuadLengthForEqualK1(arc4,0);

indQF=findcells(arc,'FamName','QF') ;
indQFM=findcells(arc,'FamName','QFM') ;

DK1_b=compK1(arc,indQF,indQFM);
LQ6_b=getcellstruct(arc,'Length',indQF)*length(indQF);
LQ8_b=getcellstruct(arc,'Length',indQFM)*length(indQFM);

DK1_a=compK1(arc1,indQF,indQFM);
LQ6_a=getcellstruct(arc1,'Length',indQF)*length(indQF);
LQ8_a=getcellstruct(arc1,'Length',indQFM)*length(indQFM);


disp('      before     after')
disp([['DK1  ';'LQF  ';'LQFM ';'sum L'],...
    num2str([[DK1_b;LQ6_b(1);LQ8_b(1);LQ6_b(1)+LQ8_b(1)],...
    [DK1_a;LQ6_a(1);LQ8_a(1);LQ6_a(1)+LQ8_a(1)]],'   %1.5f')])

