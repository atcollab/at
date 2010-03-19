function [nuh,nuv,X,Y,orh,orv]=plotnonlinmap(machine,X,Y,nbturn,type,DP,fname)
%Computes and plots the frequency map for a given lattice
%INPUTS: MACHINE is the at structure representing the accelerator with no
%time dependent fields (no cavities
%X is a vector containing the initial conditions for the horizontal
%position
%Y is a vector containing the initial position for the vertical position
%the number of initial considions tested is size(X)*size(Y). for each X all
%different Y will be tested (ans vis-versa)
%(if type=1) or the momentum deviation (if type=2).
%NBTURN is the number of turns for the traking.Traking takes place for
%2*NBTURN and each half is processed separately.
%TYPE should be 1 for an XY map, 2 for a XDP map.
%DP relevant only if TYPE=1. Defines the momentum deviation for which the
%XY map is calculated.
%FNAME Filename in which the results will be saved. If not enteed, the
%results are not saved.

%OUTPUTS:NUH and NUV are matrices of sizes: size(x)*size(Y) 4
% each lines corresponds to a given initial conditions.
% The initial conditions are aligned in the following way:(X1,Y1)
% (X2,Y1)...(Xn,Y1) (X1,Y2) (X2,Y2)...(Xn,Y2).....(X1,Yn)...(Xn,Yn).
%the first row is the tune computed on the first NBTURNS turns.
%the second row is te tune for the next NBTURN turns.
%the third row is the difference of the first two rows
%the fourth row is the DC orbit.
%H refers to horizontal, V to vertical.
% The initial conditions are aligned in the following way:(X1,Y1)
% (X2,Y1)...(Xn,Y1) (X1,Y2) (X2,Y2)...(Xn,Y2).....(X1,Yn)...(Xn,Yn).M
tic
%computes the frequency map
[nuh,nuv,X,Y,orh,orv,nu,ksi]=nonlinmap(machine,X,Y,nbturn,type,DP);
%plots the results
plotmap(X,Y,nuh,nuv);
%Saves the result in a matlab variable
if nargin==7
    
    save(fname,'nuh','nuv','X','Y','nu','ksi');
end
   toc 