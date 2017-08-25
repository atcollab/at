function [tune, varargout] = tunechrom(RING,DP,varargin)
%TUNECHROM computes linear tunes and chromaticities for COUPLED or UNCOUPLED lattice
%
% TUNE = TUNECHROM(RING,DP) - quick calculation of fractional part of the tune 
%    from numerically computed transfer matrix, assuming NO X-Y coupling.
%    If the tune is above half-integer TUNECHROM finds 1/2 - nu 
%
% TUNE = TUNECHROM(RING,DP,TUNEGUESS) - resolves the integer and half-integer 
%    uncertainty using the TUNEGUESS value. TUNEGUESS = [NUX,NUY]
%   
% [TUNE, CHROM] = TUNECHROM(RINGD,DP,TUNEGUESS,'chrom',DDP) - optionally computes
%    chromaticity by numerical differentiation from the difference between tune
%    values at momentums DP+DDP and DP
% 
% [TUNE, CHROM] = TUNECHROM(RINGD,DP,TUNEGUESS,'chrom') same as above, only uses
%    for DDP the value set in global structure NUMDIFPARAMS.
%    If NUMDIFPARAMS is not defined, TUNECHROM uses the internal default value for DDP (1e-8). 
%
% TUNECHROM(..., 'coupling') - when 'coupling' switch is added to any of the above 
%    syntax options, the tunes and chromaticities are calculated assuming
%    COUPLED lattice for two transverse eigenmodes.
%
% Note: TUNECHROM computes tunes and chromaticities from the 4-by-4
%   transfer matrix. The transfer matrix is found in FINDM44 using
%   numerical differentiation. The error of numerical differentiation 
%   is sensitive to the step size. (Reference: Numerical Recipes)
%   Calculation of tunes in TUNECHROM involves one numerical differentiation
%   to find the 4-by-4 transfer matrix.
%   Calculation of chromaticity in TUNECHROM involves TWO!!! numerical differentiations.
%   The error in calculated chromaticity from may be substantial (~ 1e-5).
%   Use the DDP argument to control the step size in chromaticity calculations
%   Another  way to control the step size is NUMDIFPARAMS structure
%   
%   
% See also LINOPT, TWISSRING, TWISSLINE, NUMDIFPARAMS

DDP_default = 1e-8;

% Process input arguments
if nargin>2
    % See if 'coupling' switch is thrown as the last argument
    if ischar(varargin{end}) & strncmp(lower(varargin{end}),'coupl',5)
        COUPLINGFLAG  = 1;
    else
        COUPLINGFLAG  = 0;
    end
    % See if TUNEGUESS is specified as the third argument
    if isnumeric(varargin{1}) & length(varargin{1})==2
        TUNEGUESSFLAG = 1;
        TUNEGUESS = varargin{1};
    else
        TUNEGUESSFLAG = 0;
        TUNEGUESS = [0.25, 0.25]; % if no TUNEGUESS is specified
    end
    % See if any of the argument is 'chrom' ,then chech if the argument after 'chrom' is DDP
    CHROMFLAG = 0;
    for i = 1:nargin-2
       if strcmp(lower(varargin{i}),'chrom')
            CHROMFLAG = 1;
            if i<nargin-2 & isnumeric(varargin{i+1})
                DDP = varargin{i+1};
            else
                % Check if NUMDIFPARAMS is defined globally
                global NUMDIFPARAMS
                if isfield(NUMDIFPARAMS,'DPStep')
                    DDP = NUMDIFPARAMS.DPStep;
                else % use default DDP
                    DDP =  DDP_default; 
                end           
            end
            break
        end
    end
            
    
else
    COUPLINGFLAG = 0;
    CHROMFLAG = 0;
    TUNEGUESSFLAG = 0;
    TUNEGUESS = [0.25, 0.25]; % if no TUNEGUESS is specified
end

M44 = findm44(RING,DP);

if COUPLINGFLAG
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
    
    cos_mu_x = trace(A)/2;
    cos_mu_y = trace(B)/2;
    
    sin_mu_x = sign(A(1,2))*sqrt(-A(1,2)*A(2,1)-(A(1,1)-A(2,2))^2/4);
    sin_mu_y = sign(B(1,2))*sqrt(-B(1,2)*B(2,1)-(B(1,1)-B(2,2))^2/4);
    
else
    cos_mu_x = (M44(1,1)+M44(2,2))/2;
    cos_mu_y = (M44(3,3)+M44(4,4))/2;    
end
TUNE = acos([cos_mu_x,cos_mu_y])/2/pi;


if TUNEGUESSFLAG
    % Check if the TUNE is in the same quadrant as TUNEGUESS
    guess_quadrant = (TUNEGUESS-floor(TUNEGUESS))> 1/2;
    tune = floor(TUNEGUESS) +  guess_quadrant + TUNE.*(sign(1/2 - guess_quadrant));
else
    tune = TUNE;
end

if CHROMFLAG & nargout > 1  
    if COUPLINGFLAG
        tune_DDP = tunechrom(RING,DP+DDP,TUNEGUESS,'coupling');
    else
        tune_DDP = tunechrom(RING,DP+DDP,TUNEGUESS);
    end
    varargout{1} = (tune_DDP - tune)/DDP;
end