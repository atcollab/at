function rmatch=MatchLast2CorForFirstBPM(ring,inCOD,indBPM,indHCor,indVCor)
%takes the last two correctors to match the orbit and angle trajectory at
%the first BPM.

% get trajectory
[t0]=findtrajectory6Err(ring,indBPM,inCOD);
 
% trn4=linepass([ring; ring; ring; ring],[inCOD(1:4),0,0]',[indBPM length(ring)+indBPM  length(ring)*2+indBPM length(ring)*3+indBPM]);
%  figure('name','initial matched 4 turn'); plot(trn4');ylim([-3e-3 3e-3]);
% export_fig('InitialTrajectory4turnRF.jpg')

%figure('name','initial'); plot(t0');ylim([-3e-3 3e-3]);
     
% match angle and position at BPM 1 of the un rotated lattice ( BPM n of the rotated
% lattice) to be identical to those of the initial trajectory
h1=atVariableBuilder(ring,indHCor(end-1),{'PolynomB',{1,1}});
h2=atVariableBuilder(ring,indHCor(end),{'PolynomB',{1,1}});
v1=atVariableBuilder(ring,indVCor(end-1),{'PolynomA',{1,1}});
v2=atVariableBuilder(ring,indVCor(end),{'PolynomA',{1,1}});
Variab=[h1 h2 v1 v2];

bpmmatchind=1;

Constr=struct(...
    'Fun',@(r,~,~)transpose(findtrajectory6Err([r;r],length(r)+indBPM(bpmmatchind),inCOD)),... % bpmmatchind BPM of second turn
    'Weight',ones(6,1)',...
    'RefPoints',[1],...
    'Min',t0(:,bpmmatchind)',...
    'Max',t0(:,bpmmatchind)');

% input optics and COD
[intwi,~,~]=atlinopt(ring,0,1);
intwi.ClosedOrbit=inCOD(1:4);

[rmatch]=atmatch(...
     ring,Variab,Constr,1e-12,100,0,@lsqnonlin,intwi);
 
 [tm]=findtrajectory6Err(rmatch,indBPM,inCOD);
%figure('name','rotated matched'); plot(tm');ylim([-3e-3 3e-3]);


% % check multi turns.
% ring=rmatch;
% trn4=linepass([ring; ring; ring; ring],[inCOD(1:4),0,0]',[indBPM length(ring)+indBPM  length(ring)*2+indBPM length(ring)*3+indBPM]);
%  figure('name','initial matched 4 turn'); plot(trn4');ylim([-3e-3 3e-3]);
% export_fig('MatchedTrajectory4turnRF.jpg')

% % check that now a COD exists.
% ring=rmatch;
%findorbit4Err(ring,0,indBPM,[inCOD 0 0]');

% % 6D ok if cavity is correctly set before setting errors in the lattice
%findorbit6Err(atsetRFCavity(ring,6.5e6,0,992,0.0),indBPM,[inCOD 0 0]'); 

return
