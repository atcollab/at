function [ringdata,elemdata] = atlinopt4(ring,varargin)
%ATLINOPT4 Performs the 4D linear analysis of COUPLED lattices
%
%[RINGDATA,ELEMDATA] = ATLINOPT4(RING,REFPTS)
%
% IMPORTANT!!! ATLINOPT4 assumes a constant momentum deviation.
%   PassMethods used for any element in the RING SHOULD NOT
%   1.change the longitudinal momentum dP
%     (cavities , magnets with radiation, ...)
%   2.have any time dependence (localized impedance, fast kickers, ...)
%
%RINGDATA is a structure array with fields:
%   tune          1x2 tunes
%   chromaticity  1x2 chromaticities (only with get_chrom or get_w flags)
%
%ELEMDATA is a structure array with fields:
%   SPos        - longitudinal position [m]
%   ClosedOrbit - 4x1 closed orbit vector with components
%                 x, px, y, py (momentums, NOT angles)
%   Dispersion  - [eta_x; eta'_x; eta_y; eta'_y] 4x1 dispersion vector
%   M           - 4x4 transfer matrix M from the beginning of RING
%                 to the entrance of the element [2]
%   A           - 2x2 matrix A in [4]
%   B           - 2x2 matrix B in [4]
%   C           - 2x2 matrix C in [4]
%   gamma       - gamma parameter of the transformation to eigenmodes [4]
%   beta        - [betax, betay] vector
%   alpha       - [alphax, alphay] vector
%   mu          - [mux, muy] Betatron phase advances
%   W           - [Wx, Wy]  Chromatic amplitude function [3] (only with the
%                           get_w flag)
% 
%   Use the Matlab function "cat" to get the data from fields of ELEMDATA as MATLAB arrays.
%   Example: 
%   >> [ringdata, elemdata] = ATLINOPT4(ring,1:length(ring));
%   >> beta = cat(1,elemdata.beta);
%   >> s = cat(1,elemdata.SPos);
%   >> plot(S,beta)
%
%   All values are specified at the entrance of each element specified in REFPTS.
%   REFPTS is an array of increasing indexes that  select elements
%   from the range 1 to length(LINE)+1. Defaults to 1 (initial point)
%   See further explanation of REFPTS in the 'help' for FINDSPOS
%
% [...] = ATLINOPT4(...,'get_chrom')
%   Trigger the computation of chromaticities
%
% [...] = ATLINOPT4(...,'get_w')
%   Trigger the computation of chromatic amplitude functions (time consuming)
%
% [...] = ATLINOPT4(...,'dp',DPP)
%   Analyses the off-momentum lattice by specifying the central
%   off-momentum value
%
% [...] = ATLINOPT4(...,'ct',DCT)
%   Analyses the off-momentum lattice by specifying the path lengthening
%   corresponding to a frequency shift. The resulting deltap/p is returned
%   in the 5th component of the ClosedOrbit field.
%
% [...] = ATLINOPT4(...,'orbit',ORBITIN)
%   Do not search for closed orbit. Instead ORBITIN,a 6x1 vector
%   of initial conditions is used: [x0; px0; y0; py0; DP; 0].
%   The sixth component is ignored.
%   This syntax is useful to specify the entrance orbit if RING is not a
%   ring or to avoid recomputing the closed orbit if is already known.
%
% [...] = ATLINOPT4(...,'twiss_in',TWISSIN)
%   Computes the optics for a transfer line.
%
% TWISSIN is a scalar structure with fields:
%   ClosedOrbit - 4x1 initial closed orbit. Default: zeros(4,1)
%   Dispersion  - 4x1 initial dispersion.   Default: zeros(4,1)
%   mu          - [ mux, muy] horizontal and vertical betatron phase
%   beta        - [betax0, betay0] vector
%   alpha       - [alphax0, alphay0] vector
%
%  REFERENCES
%	[1] D.Edwards,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
%	[2] E.Courant, H.Snyder
%	[3] Brian W. Montague Report LEP Note 165, CERN, 1979
%	[4] D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams, vol.2 (1999)
%
%  See also atlinopt2, atlinopt6

% [...] = ATLINOPT4(...,'coupled',flag)
%   Private keyword for computation without coupling (used for atlinopt2)
%
% [...] = ATLINOPT4(...,'mkey',MKEY)
%   Private keyword for the generation of the legacy atlinopt output
%
%  See also atlinopt atlinopt2 atlinopt6 tunechrom

%check_radiation(ring,false);
NE = length(ring);
[get_chrom,varargs]=getflag(varargin,'get_chrom');
[get_w,varargs]=getflag(varargs,'get_w');
[twiss_in,varargs]=getoption(varargs,'twiss_in',[]);
[orbitin,varargs]=getoption(varargs,'orbit',[]);
[ct,varargs]=getoption(varargs,'ct',NaN);
[dp,varargs]=getoption(varargs,'dp',0);
[coupled,varargs]=getoption(varargs,'coupled',true);
[mkey,varargs]=getoption(varargs,'mkey','M');
[DPStep,varargs]=getoption(varargs,'DPStep');
[XYStep,varargs]=getoption(varargs,'XYStep');

if length(varargs) >= 1
    if islogical(varargs{1})
        refpts=varargs{1}(:);
        refpts(end+1:NE+1)=false;
    elseif isnumeric(varargs{1})
        refpts=setelems(false(NE+1,1),varargs{1});
    else
        error('REFPTS must be numeric or logical');
    end
else
    refpts=setelems(false(NE+1,1),1);
end

if isempty(twiss_in)        % Circular machine
    [orbit,orbitin]=findorbit(ring,refpts,'dp',dp,'ct',ct,'orbit',orbitin,'XYStep',XYStep);
    dp=orbitin(5);
    [orbitP,o1P]=findorbit4(ring,dp+0.5*DPStep,refpts,orbitin,'XYStep',XYStep);
    [orbitM,o1M]=findorbit4(ring,dp-0.5*DPStep,refpts,orbitin,'XYStep',XYStep);
else                        % Transfer line
    if ~isempty(orbitin), orbitin=zeros(6,1); end
    orbit=linepass(ring,orbitin,refpts);
    try
        disp0=twiss_in.Dispersion;
    catch
        disp0=zeros(4,1);
    end
    dorbit=0.5*[DPStep*disp0;DPStep;0];
    orbitP=linepass(ring,orbitin+dorbit,refpts);
    orbitM=linepass(ring,orbitin-dorbit,refpts,'KeepLattice');
end
dispersion = (orbitP-orbitM)/DPStep;

[M44, MS] = findm44(ring,dp,refpts,'orbit',orbitin,'XYStep',XYStep);
T12=squeeze(num2cell(MS,[1 2]));

% Calculate A,B,C, gamma at the first element
M =M44(1:2,1:2);
N =M44(3:4,3:4);

if coupled
    m =M44(1:2,3:4);
    n =M44(3:4,1:2);
    
    % 2-by-2 symplectic matrix
    S = [0 1; -1 0];
    H = m + S*n'*S';
    t = trace(M-N);
    
    g2=(1 + sqrt(t*t/(t*t+4*det(H))))/2;
    g = sqrt(g2);
    C = -H*sign(t)/(g*sqrt(t*t+4*det(H)));
    A = g2*M  -  g*(m*S*C'*S' + C*n) + C*N*S*C'*S';
    B = g2*N  +  g*(S*C'*S'*m + n*C) + S*C'*S'*M*C;
    [MSA,MSB,gamma,CL,AL,BL]=cellfun(@analyze,T12,'UniformOutput',false);
    MSA=cat(3,MSA{:});
    MSB=cat(3,MSB{:});
    %   [BX,AX,~]=cellfun(@closure,AL);
    %   [BY,AY,~]=cellfun(@closure,BL);
else
    A = M;
    B = N;
    MSA=MS(1:2,1:2,:);
    MSB=MS(3:4,3:4,:);
    gamma=1;
    CL=zeros(2,2);
    AL=NaN(2,2);
    BL=NaN(2,2);
end

if isempty(twiss_in)            % Circular machine
    [beta0_a,alpha0_a,tune_a]=closure(A);
    [beta0_b,alpha0_b,tune_b]=closure(B);
else                            % Transfer line
    beta0_a=twiss_in.beta(1);
    beta0_b=twiss_in.beta(2);
    alpha0_a=twiss_in.alpha(1);
    alpha0_b=twiss_in.alpha(2);
    tune_a=NaN;
    tune_b=NaN;
end

[BX,AX,MX]=lop(reshape(MSA,2,2,[]),beta0_a,alpha0_a);
[BY,AY,MY]=lop(reshape(MSB,2,2,[]),beta0_b,alpha0_b);

%tunes = [MX(end),MY(end)]/2/pi;

ringdata=struct('tune',{[tune_a tune_b]});

elemdata = struct(...
    'SPos',num2cell(findspos(ring,refpts))',...
    'ClosedOrbit',num2cell(orbit(1:4,:),1)',...
    'Dispersion',num2cell(dispersion(1:4,:),1)',...
    mkey,T12,...
    'gamma',gamma, 'C',CL, 'A',AL, 'B',BL,...
    'beta',num2cell([BX,BY],2),...
    'alpha',num2cell([AX,AY],2),...
    'mu',num2cell([MX,MY],2));

if get_w
    % Calculate optics for DP +/- DDP
    [rp,lp] = atlinopt4(ring,refpts,'dp',dp+0.5*DPStep,'orbit',o1P,'coupled',coupled);
    [rm,lm] = atlinopt4(ring,refpts,'dp',dp-0.5*DPStep,'orbit',o1M,'coupled',coupled);
    w = chromfunc(DPStep,cat(1,lp.alpha),cat(1,lp.beta),cat(1,lm.alpha),cat(1,lm.beta));
    [elemdata.W]=deal(w{:});
    ringdata.chromaticity=(rp.tune - rm.tune)/DPStep;
elseif get_chrom
    % Calculate tunes for DP +/- DDP
    tunep=tunechrom(ring,dp+0.5*DPStep,'orbit',o1P,'coupled',coupled);
    tunem=tunechrom(ring,dp-0.5*DPStep,'orbit',o1M,'coupled',coupled);
    ringdata.chromaticity=(tunep-tunem)/DPStep;
end

    function [E12,F12,gamma,CL,AL,BL]=analyze(MS)
        M12 =MS(1:2,1:2);
        N12 =MS(3:4,3:4);
        m12 =MS(1:2,3:4);
        n12 =MS(3:4,1:2);
        
        gamma = sqrt(det(n12*C+g*N12));
        E12 = (g*M12-m12*S*C'*S')/gamma;
        F12 = (n12*C+g*N12)/gamma;
        
        CL=(M12*C+g*m12)*S*F12'*S';
        AL=E12*A*S*E12'*S';
        BL=F12*B*S*F12'*S';
    end

    function mask=setelems(mask,idx)
        mask(idx)=true;
    end

    function UP = BetatronPhaseUnwrap(P)
        JUMPS = diff([0;P],1,1) < -1.e-3;
        UP = P+cumsum(JUMPS)*2*pi;
    end

    function [beta,alpha,phase]=lop(MS,beta0,alpha0)
        bbb=squeeze(MS(1,2,:));
        aaa=squeeze(MS(1,1,:))*beta0-bbb*alpha0;
        
        beta = (aaa.^2 + bbb.^2)/beta0;
        alpha = -(aaa.*squeeze(MS(2,1,:)*beta0-MS(2,2,:)*alpha0) + bbb.*squeeze(MS(2,2,:)))/beta0;
        try
            phase = atan2(bbb,aaa);
        catch
            phase=NaN(size(beta));
        end
        phase = BetatronPhaseUnwrap(phase);
    end

    function [beta,alpha,tune] = closure(AB)
        cosmu = (AB(1,1) + AB(2,2))/2;
        diff  = (AB(1,1) - AB(2,2))/2;
        sinmu = sign(AB(1,2))*sqrt(-AB(1,2)*AB(2,1)-diff*diff);
        try
            alpha = (AB(1,1)-AB(2,2))/2/sinmu;
            beta = AB(1,2)/sinmu;
            tune = mod(atan2(sinmu,cosmu)/2/pi,1);
        catch           % Unstable ring
            alpha = NaN;
            beta = NaN;
            tune = NaN;
        end
    end

    function w=chromfunc(ddp,aup,bup,adn,bdn)
        % Compute the chromatic W function
        db = (bup - bdn) / ddp;
        mb = (bup + bdn) / 2;
        da = (aup - adn) / ddp;
        ma = (aup + adn) / 2;
        w = num2cell(sqrt((da - ma ./ mb .* db).^2 + (db ./ mb).^2),2);
    end

end
