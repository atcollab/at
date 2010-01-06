function [ENVELOPE, RMSDP, RMSBL] = ohmienvelope(RING,RADELEMINDEX,varargin)
%OHMIENVELOPE calculates equilibrium beam envelope in a 
% circular accelerator using Ohmi's beam envelope formalism [1]
% [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
% 
% [ENVELOPE, RMSDP, RMSBL] = ONMIENVELOPE(RING,RADELEMINDEX,REFPTS)
% 
% ENVELOPE is a structure with fields
% Sigma   - [SIGMA(1); SIGMA(2)] - RMS size [m] along 
%           the principal axis of a tilted ellips 
%           Assuming normal distribution exp(-(Z^2)/(2*SIGMA))
% Tilt    - Tilt angle of the XY ellips [rad]
%           Positive Tilt corresponds to Corkscrew (right) 
%           rotatiom of XY plane around s-axis
% R       - 6-by-6 equilibrium envelope matrix R
%
% RMSDP   - RMS momentum spread 
% RMSBL   - RMS bunch length[m]

NumElements = length(RING);
 
[MRING, MS, orbit] = findm66(RING,1:NumElements+1);

B = cell(1,NumElements); % B{i} is the diffusion matrix of the i-th element
BATEXIT = B;             % BATEXIT{i} is the cumulative diffusion matrix from 
                         % element 0 to the end of the i-th element
M = B;                   % 6-by-6 transfer matrix of the i-th element

% calculate Radiation-Diffusion matrix B for elements with radiation
for i = RADELEMINDEX
   B{i} = findmpoleraddiffmatrix(RING{i},orbit(:,i));
end

% Calculate 6-by-6 linear transfer matrix in each element
% near the equilibrium orbit 
for i = 1:NumElements
   ELEM = RING{i};
   M{i} = findelemm66(ELEM,ELEM.PassMethod,orbit(:,i));
   % Set Radiation-Diffusion matrix B to 0 in elements without radiation
   if isempty(B{i})
      B{i} = zeros(6,6);
   end
end
% Calculate cumulative Radiation-Diffusion matrix for the ring
BCUM = zeros(6,6); % Cumulative diffusion matrix of the entire ring

for i = 1:NumElements
   BCUM = M{i}*BCUM*M{i}' + B{i};
   BATEXIT{i} = BCUM;
end
% ------------------------------------------------------------------------
% Equation for the moment matrix R is
%         R = MRING*R*MRING' + BCUM;
% We rewrite it in the form of Lyapunov equation to use MATLAB's LYAP function
%            AA*R + R*BB = -CC  
% where 
%				AA =  inv(MRING)
%				BB = -MRING'
%				CC = -inv(MRING)*BCUM
% -----------------------------------------------------------------------
AA =  inv(MRING);
BB = -MRING';
CC = -inv(MRING)*BCUM;
 
R = lyap(AA,BB,CC);     % Envelope matrix at the rinng entrance

RMSDP = sqrt(R(5,5));   % R.M.S. energy spread
RMSBL = sqrt(R(6,6));   % R.M.S. bunch lenght

if nargin == 2 % no REFPTS
    ENVELOPE.R = R;
elseif nargin == 3
    REFPTS = varargin{1};
    
    REFPTSX = REFPTS;
    % complete REFPTS with 1 and NumElements+1 if necessary
    if REFPTS(1)~=1
        REFPTSX = [1 REFPTS];
    end
    if REFPTS(end)~= NumElements+1
        REFPTSX = [REFPTSX NumElements+1];
    end
    % Now REFPTS has at least 2 ponts and the first one is the ring entrance
    
    NRX = length(REFPTSX);
    ENVELOPE = struct('Sigma',num2cell(zeros(2,NRX),1),...
        'Tilt',0,'R',zeros(6)); 
    
    ENVELOPE(1).R = R;

    for i=2:NRX
        ELEM = REFPTSX(i);
        ENVELOPE(i).R = MS(:,:,ELEM)*R*MS(:,:,ELEM)'+BATEXIT{ELEM-1};
    end
    
   
    if REFPTS(1)~=1
        ENVELOPE = ENVELOPE(2:end);
    end
    if REFPTS(end)~= NumElements+1
        ENVELOPE = ENVELOPE(1:end-1);
    end

else
    error('Too many input arguments');
end

for i=1:length(ENVELOPE)
    [U,DR] = eig(ENVELOPE(i).R([1 3],[1 3]));
    ENVELOPE(i).Tilt = asin((U(2,1)-U(1,2))/2);
    ENVELOPE(i).Sigma(1) = sqrt(DR(1,1));
    ENVELOPE(i).Sigma(2) = sqrt(DR(2,2));
end