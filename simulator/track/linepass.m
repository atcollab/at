function Rout = linepass(LINE,Rin,varargin)
%LINEPASS tracks particles through a sequence of elements 
% in the cell array LINE. For each element it calls 
% the pass-method specified in the 'PassMethod' field. 
%		
% Rout = LINEPASS(LINE,Rin) tracks particle(s) to the end 
%     of the LINE. Rin and Rout are 6-by-1 
%     column vectors or 6-by-N matrixes
%     Each column represents a different initial condition or particle.
%
% Rout = LINEPASS(LINE,Rin,REFPTS) also returns intermediate results 
%     at the entrance of each element specified in the REFPTS
%
%     REFPTS is an array of increasing indexes that  selects elements 
%     between 1 and length(LINE)+1. 
%     See further explanation of REFPTS in the 'help' for FINDSPOS
%     
%     NOTE:
%     LINEPASS(LINE,Rin,length(LINE)+1) is the same as  LINEPASS(LINE,Rin)
%     since the reference point length(LINE)+1 is the exit of the last element
%     LINEPASS(LINE,Rin,1) is a copy of Rin since the 
%     reference point 1 is the entrance of the first element
%     
%     OUTPUT FORMAT:
%     Rout is 6-by-(number of columns in Rin)*length(REFPTS) matrix
%     where blocks 6-by-(number of columns in Rin) corresponds 
%     to different REFPTS
%     FOR EXAMPLE:
%     if Rin is 6-by-2 maid of two 6-by-1 column vectors [Rin1, Rin2]
%     and REFPTS = [N1 N2 N3] so that N1<N2<N3
%     the output is [Rout1(N1) Rout2(N1) Rout1(N2) Rout2(N2) Rout1(N3) Rout2(N3)]
%
% LINEPASS(LINE,Rin,REFPTS,'reuse') and  LINEPASS(LINE,Rin,'reuse')
%    with 'reuse' flag is more efficient because
%    it reuses some of the data and functions stored in 
%    persistent memory from the previous calls to LINEPASS. 
%
%    !!! In order to use this option, LINEPASS must first be
%    called without the reuse flag. This will 
%    create persistent data structures and keep pointers 
%    to pass-method functions. 
%
%    !!! LINEPASS(...'reuse') assumes that the number of 
%    elements in LINE and pass methods specified in the
%    PassMethod field of each element DO NOT CHANGE between
%    calls. Otherwise, LINEPASS without 'reuse' must 
%    be called again. The values of elements fields such as 'Length' or
%    'K' are allowed to change
%     
% See also: RINGPASS FINDSPOS


test = strcmpi(varargin,'reuse');
if any(test)
    NEWLATTICEFLAG = 0;
else 
    NEWLATTICEFLAG = 1;
end

numericargs = varargin(find(~test));


if length(numericargs) > 0
    REFPTS = numericargs{1};
else
    REFPTS = length(LINE)+1;
end 
    
Rout = atpass(LINE,Rin,NEWLATTICEFLAG,1,REFPTS);

    