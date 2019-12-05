function [TD, varargout] = twissring(RING,DP,varargin)
%TWISSRING calculates linear optics functions for an UNCOUPLED ring
% 
% [TwissData, tune]  = TWISSRING(LATTICE,DP) calculates twiss parameters
%    and closed orbit coordinates at the RING entrance assuming
%    constant energy deviation DP.
%
% [TwissData, tune]  = TWISSRING(LATTICE,DP,REFPTS) calculates Twiss parameters
%    and closed orbit coordinates at specified reference points REFPTS.
%
%    Note: REFPTS is an array of increasing indexes that  
%    select elements from range 1 to length(LATTICE)+1. 
%    See further explanation of REFPTS in the 'help' for FINDSPOS 
%
% [TwissData, tune, chrom]  = TWISSRING(...,'chrom', DDP) also calculates
%    linear dispersion and chromaticity. Dispersion is returned as one 
%    of the fields in TwissData.
%    !!! Last argument DDP is a momentum deviation on top
%    of DP (the second argument) used to calculate and normalize
%    dispersion and chromaticity. If not supplied
%    the default value of 1e-8 is used.
%
%    Note: To resolve the integer part of the tune
%    and the uncertainty of acos(trace(M)/2) it is necessary to
%    supply sufficient number of REFPTS properly spaced in betatron phase.
%
% TwisData is a 1-by-REFPTS (1-by-1) structure array with fields
%       (Some are the same as in the output of LINOPT)
%       ElemIndex   - integer (element number) in the RING 
%       SPos        - longitudinal position [m]
%       ClosedOrbit - closed orbit column vector with 
%                     components x, px, y, py (momentums, NOT angles)						
%       Dispersion  - dispersion orbit position 4-by-1 vector with 
%                     components [eta_x, eta_prime_x, eta_y, eta_prime_y]'
%                     calculated with respect to the closed orbit with 
%                     momentum deviation DP
%       M44         - 4x4 transfer matrix M from the beginning of RING
%                     to the entrance of the element for specified DP [2]
%       beta        - [betax, betay] horizontal and vertical Twiss parameter beta
%       alpha       - [alphax, alphay] horizontal and vertical Twiss parameter alpha
%       mu          - [mux, muy] horizontal and vertical betatron phase
%                     !!! NOT 2*PI normalized
% 
% Use MATLAB function CAT to get the data from fields of TwissData into MATLAB arrays.
%     Example: 
%     >> TD = twissring(THERING,0,1:length(THERING));
%     >> BETA = cat(1,TD.beta);
%     >> S = cat(1,TD.SPos);
%     >> plot(S,BETA(:,1))
%  
% See also TWISSLINE, LINOPT, TUNECHROM.

NE=length(RING);
% Process input arguments
[CHROMFLAG,args]=getflag(varargin,'chrom');
[REFPTS,DDP]=getargs(args,{NE+1,1.e-8});
CHROMFLAG=CHROMFLAG || (nargout == 3);

% Include the endpoint if it is not already in REFPTS
if REFPTS(end)==NE+1
    [M44, MS, orb] = findm44(RING,DP,REFPTS);
else
    [M44, MS, orb] = findm44(RING,DP,[REFPTS,NE+1]);
end

cos_mu_x = (M44(1,1)+M44(2,2))/2;
cos_mu_y = (M44(3,3)+M44(4,4))/2;

sin_mu_x = sign(M44(1,2))*sqrt(-M44(1,2)*M44(2,1)-(M44(1,1)-M44(2,2))^2/4);
sin_mu_y = sign(M44(3,4))*sqrt(-M44(3,4)*M44(4,3)-(M44(3,3)-M44(4,4))^2/4);


ax = (M44(1,1)-M44(2,2))/2/sin_mu_x;
ay = (M44(3,3)-M44(4,4))/2/sin_mu_y;

bx = M44(1,2)/sin_mu_x;
by = M44(3,4)/sin_mu_y;

BX = squeeze((MS(1,1,:)*bx-MS(1,2,:)*ax).^2 + MS(1,2,:).^2)/bx;
BY = squeeze((MS(3,3,:)*by-MS(3,4,:)*ay).^2 + MS(3,4,:).^2)/by;


AX = -squeeze((MS(1,1,:)*bx-MS(1,2,:)*ax).*(MS(2,1,:)*bx-MS(2,2,:)*ax) + MS(1,2,:).*MS(2,2,:))/bx;
AY = -squeeze((MS(3,3,:)*by-MS(3,4,:)*ay).*(MS(4,3,:)*by-MS(4,4,:)*ay) + MS(3,4,:).*MS(4,4,:))/by;

MX = atan2(squeeze(MS(1,2,:)), squeeze(MS(1,1,:)*bx-MS(1,2,:)*ax));
MY = atan2(squeeze(MS(3,4,:)), squeeze(MS(3,3,:)*by-MS(3,4,:)*ay));

MX = BetatronPhaseUnwrap(MX);
MY = BetatronPhaseUnwrap(MY);

tune=mod(atan2([sin_mu_x,sin_mu_y],[cos_mu_x cos_mu_y])/2/pi,1.0);

NR = length(REFPTS);
% Build TD only for points originally referenced in REFPTS
TD = struct('ElemIndex',num2cell(REFPTS),...
    'SPos',num2cell(findspos(RING,REFPTS)),...
    'ClosedOrbit',num2cell(orb(:,1:NR),1),...
    'M44', squeeze(num2cell(MS(:,:,1:NR),[1 2]))',...
    'beta', num2cell([BX(1:NR),BY(1:NR)],2)',...
    'alpha', num2cell([AX(1:NR),AY(1:NR)],2)',...
    'mu', num2cell([MX(1:NR),MY(1:NR)],2)');


if CHROMFLAG
    [TD_DDP,tune_DDP] = twissring(RING,DP+DDP,REFPTS);
    DORBIT = reshape(cat(1,TD_DDP.ClosedOrbit),4,[]);
    DISPERSION = num2cell((DORBIT-orb(:,1:NR))/DDP,1);
    [TD.Dispersion] = deal( DISPERSION{:});
    varargout{2} = (tune_DDP-tune)/DDP;
end
    
if nargout>1
    varargout{1}=tune;
end

function UP = BetatronPhaseUnwrap(P)
    % unwrap negative jumps in betatron
    %JUMPS = [0; diff(P)] < -1.e-5;
    JUMPS = [0; diff(P)] < -1.e-3;
    UP = P+cumsum(JUMPS)*2*pi;


