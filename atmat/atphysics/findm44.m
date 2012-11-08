function [M44, varargout]  = findm44(LATTICE,DP,varargin)
%FINDM44 numerically finds the 4x4 transfer matrix of an accelerator lattice
% for a particle with relative momentum deviation DP
%
% IMPORTANT!!! FINDM44 assumes constant momentum deviation.
%   PassMethod used for any element in the LATTICE SHOULD NOT
%   1.change the longitudinal momentum dP 
%     (cavities , magnets with radiation, ...)
%   2.have any time dependence (localized impedance, fast kickers, ...) 
%
% M44 = FINDM44(LATTICE,DP) finds a full one-turn 
%    matrix at the entrance of the first element
%    !!! With this syntax FINDM44 assumes that the LATTICE 
%    is a ring and first finds the closed orbit
%    
% [M44,T] = FINDM44(LATTICE,DP,REFPTS) also returns
%    4-by-4 transfer matrixes  between entrance of 
%    the first element and each element indexed by REFPTS. 
%    T is 4-by-4-by-length(REFPTS) 3 dimensional array
%    so that the set of indexes (:,:,i) selects the 4-by-4 
%    matrix at the i-th reference point.
%    
%    Note: REFPTS is an array of increasing indexes that  
%    select elements from range 1 to length(LATTICE)+1. 
%    See further explanation of REFPTS in the 'help' for FINDSPOS 
%    When REFPTS= [ 1 2 .. ] the fist point is the entrance of the 
%    first element and T(:,:,1) - identity matrix
%    
%    Note: REFPTS is allowed to go 1 point beyond the 
%    number of elements. In this case the last point is 
%    the EXIT of the last element. If LATTICE is a RING
%    it is also the entrance of the first element 
%    after 1 turn: T(:,:,end) = M
%
% [M44, T] = FINDM44(LATTICE,DP,REFPTS,ORBITIN) - Does not search for
%   closed orbit. Instead the ORBITIN,a 1-by-6 vector of initial
%   conditions is used: [x0, px0, y0, py0, DP, 0]' where
%   the same DP as argument 2. The sixth component is ignored.
%   This syntax is useful to specify the entrance orbit
%   if LATTICE is not a ring or to avoid recomputting the
%   closed orbit if is already known.
% 
% [M44, T] = FINDM44(LATTICE,DP,REFPTS,ORBITIN,'full') - same as above except  
%    matrixes returned in T are full 1-turn matrixes at the entrance of each
%    element indexed by REFPTS. 
%
% [M44, T, orbit] = FINDM44(...) in addition returns 
%    at REFPTS the closed orbit calculated along the 
%    way with findorbit4  
%
% See also LINEPASS, FINDORBIT4 FINDSPOS

% *************************************************************************
%   The numerical differentiation in FINDM44 uses symmetric form
%
%         F(x+delta) - F(x-delta)
%       --------------------------------------
%              2*delta 
%
%    with optimal differentiation step delta given by !!!! DO LATER 
%    The relative error in the derivative computed this way 
%    is !!!!!!!!!!!!!!!!! DO LATER  
%    Reference: Numerical Recipes.   


if ~iscell(LATTICE)
   error('First argument must be a cell array'); 
end

NE = length(LATTICE);


switch nargin
case 5 % FINDM44(LATTICE,DP,REFPTS,ORBITIN,'full')
    if(lower(varargin{3})=='full')
        FULLFLAG = 1;
        REFPTS = varargin{1};
        R0 = varargin{2};
        R0(5) = DP;
        R0(6)= 0;
    else
        error('Fifth argument - unknown option')
    end
case 4 % FINDM44(LATTICE,DP,REFPTS,ORBITIN)
    FULLFLAG = 0;
    REFPTS = varargin{1};
    R0 = varargin{2};
    R0(5) = DP;
    R0(6)= 0;
case 3 % FINDM44(LATTICE,DP,REFPTS)
    FULLFLAG = 0;
   	REFPTS = varargin{1};
    R0 = [findorbit4(LATTICE,DP);DP;0];
case 2 % FINDM44(LATTICE,DP)
   	REFPTS = NE+1;
    FULLFLAG = 0;
    R0 = [findorbit4(LATTICE,DP);DP;0];
otherwise
    error('Incorrect number of input arguments');
end

NR = length(REFPTS);


% Dteremine step size to use for numerical differentiation
global NUMDIFPARAMS

if isfield(NUMDIFPARAMS,'XYStep')
    d = NUMDIFPARAMS.XYStep';
else
    % optimal differentiation step - Numerical Recipes
    d =  6.055454452393343e-006; 
end
    

% Put together matrix of initial conditions

D = d*eye(4);
% First 8 columns for derivative
% 9-th column is for closed orbit
RM = [[R0 R0 R0 R0 R0 R0 R0 R0] + [D -D; zeros(2,8)],R0];

if nargout < 2 	
   % Calculate M44 at the first element only. Use linepass
   TMAT = linepass(LATTICE,RM);
   M44 = (TMAT(1:4,1:4)-TMAT(1:4,5:8))/(2*d);
   return
else					
   % Calculate matrixes at all REFPTS. Use linepass
   % Need to include the exit of the LATTICE to REFPTS array
   if NR==0 || REFPTS(NR)~=NE+1
       NR1 = NR+1;
       REFPTS(NR1)=NE+1;
   else
       NR1 = NR;
   end

   TMAT = linepass(LATTICE,RM,REFPTS);
   TMAT3 = reshape(TMAT(1:4,:),4,9,NR1);
   M44 = (TMAT3(1:4,1:4,NR1)-TMAT3(1:4,5:8,NR1))/(2*d);
   
   MSTACK = (TMAT3(:,1:4,1:NR)-TMAT3(:,5:8,1:NR))/(2*d);
   
   if FULLFLAG
       S2 = [0 1;-1 0];
       S4 = [S2, zeros(2);zeros(2),S2]; % symplectic identity matrix
       for k =1:NR
            T =  MSTACK(:,:,k);
            varargout{1}(:,:,k) = T*M44*S4'*T'*S4;
       end
   else
       varargout{1}=MSTACK; 
   end
   % return the closed orbit if requested
   if nargout == 3
       varargout{2}=squeeze(TMAT3(:,9,1:NR));
   end
       
end

