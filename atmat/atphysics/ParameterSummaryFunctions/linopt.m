function [LinData, varargout] = linopt(RING,DP,varargin);
%LINOPT performs linear analysis of the COUPLED lattices
%   Notation is the same as in reference [3]
%
%
% LinData = LINOPT(RING,DP,REFPTS) is a MATLAB structure array with fields
%    
%   ElemIndex   - ordinal position in the RING 
%   SPos        - longitudinal position [m]
%   ClosedOrbit - closed orbit column vector with 
%                 components x, px, y, py (momentums, NOT angles)						
%   Dispersion  - dispersion orbit position vector with 
%                 components eta_x, eta_prime_x, eta_y, eta_prime_y
%                 calculated with respect to the closed orbit with 
%                 momentum deviation DP
%   M44         - 4x4 transfer matrix M from the beginning of RING
%                 to the entrance of the element for specified DP [2]
%   A           - 2x2 matrix A in [3]
%   B           - 2x2 matrix B in [3]
%   C           - 2x2 matrix C in [3]			
%   gamma       - gamma parameter of the transformation to eigenmodes 
%   mu          - [ mux, muy] horizontal and vertical betatron phase
%   beta        - [betax, betay] vector
%
%   All values are specified at the entrance of each element specified in REFPTS. 
%   REFPTS is an array of increasing indexes that  select elements 
%   from the range 1 to length(LINE)+1. 
%   See further explanation of REFPTS in the 'help' for FINDSPOS 
%
% [LinData,NU] = LINOPT() returns a vector of linear tunes
%   [nu_u , nu_v] for two normal modes of linear motion [1] 
%
% [LinData,NU, KSI] = LINOPT() returns a vector of chromaticities ksi = d(nu)/(dP/P)
%   [ksi_u , ksi_v] - derivatives of [nu_u , nu_v] 
%
% See also FINDSPOS TWISSRING TUNECHROM
%
%   [1] D.Edwars,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
%   [2] E.Courant, H.Snyder
%   [3] D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams, vol.2 (1999)

 

NE=length(RING);
if(nargin==2)
   REFPTS= 1;
else
   REFPTS=varargin{1};
end

NR=length(REFPTS);
  

spos = findspos(RING,REFPTS);
[M44, MS, orb] = findm44(RING,DP,REFPTS);

LinData = struct('ElemIndex',num2cell(REFPTS),'SPos',num2cell(spos),...
    'ClosedOrbit',num2cell(orb,1),'M44',squeeze(num2cell(MS,[1 2]))');

% Calculate A,B,C, gamma at the first element
M =M44(1:2,1:2);
N =M44(3:4,3:4);
m =M44(1:2,3:4);
n =M44(3:4,1:2);

% 2-by-2 symplectic matrix
S = [0 1; -1 0];
H = m + S*n'*S';
t = trace(M-N);

g = sqrt(1 + sqrt(t*t/(t*t+4*det(H))))/sqrt(2);
G = diag([g g]);
C = -H*sign(t)/(g*sqrt(t*t+4*det(H)));
A = G*G*M  -  G*(m*S*C'*S' + C*n) + C*N*S*C'*S';
B = G*G*N  +  G*(S*C'*S'*m + n*C) + S*C'*S'*M*C;


   
if REFPTS(1)==1 & NR>1
    START = 2;
    LinData(1).A=A;
    LinData(1).B=B;
    LinData(1).C=C;
    LinData(1).gamma=g;
    LinData(1).beta(1) = A(1,2)/sin(acos(trace(A/2)));
    LinData(1).beta(2) = B(1,2)/sin(acos(trace(B/2)));
else
    START = 1;
end


      
    


% find  matrixes in all elements indexed by REFPTS
for i=START:NR;
    M12 =LinData(i).M44(1:2,1:2);
    N12 =LinData(i).M44(3:4,3:4);
    m12 =LinData(i).M44(1:2,3:4);
    n12 =LinData(i).M44(3:4,1:2);
   
    g2 = sqrt(det(n12*C+G*N12));
    E12 = (G*M12-m12*S*C'*S')/g2;
    F12 = (n12*C+G*N12)/g2;
   
    LinData(i).gamma=g2;
    LinData(i).C=(M12*C+G*m12)*S*F12'*S';
    LinData(i).A=E12*A*S*E12'*S';
    LinData(i).B=F12*B*S*F12'*S';
    LinData(i).beta(1) = LinData(i).A(1,2)/sin(acos(trace(A/2)));
    LinData(i).beta(2) = LinData(i).B(1,2)/sin(acos(trace(B/2)));
 
end



if nargout > 1 
   cos_mu_x = trace(A)/2;
   cos_mu_y = trace(B)/2;
   varargout{1} = acos([cos_mu_x cos_mu_y])/2/pi;
end

if nargout == 3
    global NUMDIFPARAMS

    if isfield(NUMDIFPARAMS,'DPStep')
        dDP = NUMDIFPARAMS.DPStep';
    else
        dDP =  1e-8;
    end
    % Calculate tunes for DP+dDP
    [LD, TUNES] = linopt(RING,DP+dDP,1);
    varargout{2} = (TUNES - varargout{1})/dDP;
end