function [t] = findtrajectory6Err( ring,indBPM,inCOD )
%[t    6xNbpm array of  trajectory   
% ] = findtrajectory6Err( 
% ring,   1) AT lattice
% indBPM, 2) bpm indexes
% inCOD ) 3) 6x1 input coordinates (default: zero)
%
% uses linepass to obtain trajectory in ring at indBPM
% if present, BPM errors are included on x(1st) and y(3rd) coordinates.
% 
%see also: linepass bpm_matrices bpm_process

if nargin<3
    inCOD=[0 0 0 0 0 0]';
end

% linepass
outtr=linepass(ring,inCOD,indBPM);
ox=outtr(1,:);
oy=outtr(3,:);

% set bpm errors
[rel,tel,trand] = bpm_matrices(ring(indBPM));
bpmreading = bpm_process([ox; oy],rel,tel,trand);
t=outtr;
t(1,:)=bpmreading(1,:);
t(3,:)=bpmreading(2,:);

end

