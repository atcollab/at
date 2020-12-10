function [lindata, varargout] = atlinopt(RING,dp,varargin)
%ATLINOPT Performs linear analysis of the COUPLED lattices
%
% LinData = ATLINOPT(RING,DP,REFPTS) is a MATLAB structure array with fields
%
%   ElemIndex   - ordinal position in the RING
%   SPos        - longitudinal position [m]
%   ClosedOrbit - 4x1 closed orbit vector with components
%                 x, px, y, py (momentums, NOT angles)
%   Dispersion  - 4x1 dispersion with components
%                 eta_x, eta_prime_x, eta_y, eta_prime_y
%   M44         - 4x4 transfer matrix M from the beginning of RING
%                 to the entrance of the element [2]
%   A           - 2x2 matrix A in [3]
%   B           - 2x2 matrix B in [3]
%   C           - 2x2 matrix C in [3]
%   gamma       - gamma parameter of the transformation to eigenmodes
%   mu          - [mux, muy] horizontal and vertical betatron phase advances
%   beta        - [betax, betay] vector
%   alpha       - [alphax, alphay] vector
%
%   All values are specified at the entrance of each element specified in REFPTS.
%   REFPTS is an array of increasing indexes that  select elements
%   from the range 1 to length(LINE)+1. Defaults to 1 (initial point)
%   See further explanation of REFPTS in the 'help' for FINDSPOS
%
% [LinData,NU] = ATLINOPT() returns a vector of linear tunes
%   [nu_u , nu_v] for two normal modes of linear motion [1]
%
% [LinData,NU, KSI] = ATLINOPT() returns a vector of chromaticities ksi = d(nu)/(dP/P)
%   [ksi_u , ksi_v] - derivatives of [nu_u , nu_v]
%
% [...] = ATLINOPT(...,'orbit',ORBITIN)
% [...] = ATLINOPT(RING,DP,REFPTS,ORBITIN)  (Deprecated syntax)
%   Do not search for closed orbit. Instead ORBITIN,a 6x1 vector
%   of initial conditions is used: [x0; px0; y0; py0; DP; 0].
%   The sixth component is ignored.
%   This syntax is useful to specify the entrance orbit if RING is not a
%   ring or to avoid recomputing the closed orbit if is already known.
%
% [...] = ATLINOPT(...,'coupled',flag)
%   If flag is false, a faster computation is performed
%   assuming no coupling in the lattice. Default: true
%
% [...] = ATLINOPT(...,'ct',CT)
%   Instead of computing the linear optics for  the specified DP/P,
%   computes for the path lenghening specified by CT.
%   The DP argument is ignored.
%
% [...] = ATLINOPT(...,'twiss_in',TWISSIN)
%   Computes the optics for a transfer line.
%
% TWISSIN is a scalar structure with fields:
%   ClosedOrbit - 4x1 initial closed orbit. Default: zeros(4,1)
%   Dispersion  - 4x1 initial dispersion.   Default: zeros(4,1)
%   mu          - [ mux, muy] horizontal and vertical betatron phase
%   beta        - [betax0, betay0] vector
%   alpha       - [alphax0, alphay0] vector
%
% Difference with linopt: Fractional tunes 0<=tune<1
%			  Dispersion output
%			  Alpha output
%			  Phase advance output
%			  Option to skip closed orbit search
%  REFERENCES
%    [1] D.Edwars,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
%    [2] E.Courant, H.Snyder
%    [3] D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams, vol.2 (1999)
%
%  See also atx atmodl findspos twissring tunechrom

NE = length(RING);
[twiss_in,varargs]=getoption(varargin,'twiss_in',[]);
[orbitin,varargs]=getoption(varargs,'orbit',[]);
[ct,varargs]=getoption(varargs,'ct',NaN);
[coupled,varargs]=getoption(varargs,'coupled',true);
[DPStep,varargs]=getoption(varargs,'DPStep');
[XYStep,varargs]=getoption(varargs,'XYStep');

if length(varargs) >= 2
    orbitin=varargs{2};
end
if length(varargs) >= 1
    if islogical(varargs{1})
        REFPTS=varargs{1}(:);
        REFPTS(end+1:NE+1)=false;
    elseif isnumeric(varargs{1})
        REFPTS=setelems(false(NE+1,1),varargs{1});
    else
        error('REFPTS must be numeric or logical');
    end
else
    REFPTS=setelems(false(NE+1,1),1);
end

if isempty(twiss_in)        % Circular machine
    if ~isempty(orbitin)
        if length(orbitin) >= 5
            dp=orbitin(5);
        end
        orbitin=[orbitin(1:4);dp;0];
    elseif isnan(ct)
        [~,orbitin]=findorbit4(RING,dp,REFPTS);
    else
        [~,orbitin]=findsyncorbit(RING,ct,REFPTS);
        dp=orbitin(5);
    end
    [orbitP,o1P]=findorbit4(RING,dp+0.5*DPStep,REFPTS,orbitin);
    [orbitM,o1M]=findorbit4(RING,dp-0.5*DPStep,REFPTS,orbitin);
else                        % Transfer line
    if ~isempty(orbitin)
        orbitin=zeros(6,1);
    end
    try
        disp0=twiss_in.Dispersion;
    catch
        disp0=zeros(4,1);
    end
    dorbit=0.5*[DPStep*disp0;DPStep;0];
    orbitP=linepass(RING,orbitin+dorbit,REFPTS);
    orbitM=linepass(RING,orbitin-dorbit,REFPTS,'KeepLattice');
end
orbit=linepass(RING,orbitin,REFPTS);
dispersion = (orbitP-orbitM)/DPStep;

spos = findspos(RING,REFPTS);
[M44, MS] = findm44(RING,dp,REFPTS,'orbit',orbitin,'XYStep',XYStep);
T12=squeeze(num2cell(MS,[1 2]))';

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
    
    g = sqrt(1 + sqrt(t*t/(t*t+4*det(H))))/sqrt(2);
    G = diag([g g]);
    C = -H*sign(t)/(g*sqrt(t*t+4*det(H)));
    A = G*G*M  -  G*(m*S*C'*S' + C*n) + C*N*S*C'*S';
    B = G*G*N  +  G*(S*C'*S'*m + n*C) + S*C'*S'*M*C;
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

if nargout >= 2
    varargout{1}=[tune_a tune_b];
end

if nargout >= 3
    % Calculate tunes for DP+dDP
    [LD, tunesP] = atlinopt(RING,dp+0.5*DPStep,[],'orbit',o1P,'coupled',coupled); %#ok<ASGLU>
    [LD, tunesM] = atlinopt(RING,dp-0.5*DPStep,[],'orbit',o1M,'coupled',coupled); %#ok<ASGLU>
%     tunesP=tunechrom(RING,dp+0.5*DPStep,'orbit',o1P,'coupled',coupled);
%     tunesM=tunechrom(RING,dp-0.5*DPStep,'orbit',o1M,'coupled',coupled);
    varargout{2} = (tunesP - tunesM)/DPStep;
end

lindata = struct('ElemIndex',num2cell(find(REFPTS))',...
    'SPos',num2cell(spos),...
    'ClosedOrbit',num2cell(orbit(1:4,:),1),...
    'Dispersion',num2cell(dispersion(1:4,:),1),...
    'M44',T12,...
    'gamma',gamma, 'C',CL, 'A',AL, 'B',BL,...
    'beta',num2cell([BX,BY],2)',...
    'alpha',num2cell([AX,AY],2)',...
    'mu',num2cell([MX,MY],2)');

    function [E12,F12,gamma,CL,AL,BL]=analyze(MS)
        M12 =MS(1:2,1:2);
        N12 =MS(3:4,3:4);
        m12 =MS(1:2,3:4);
        n12 =MS(3:4,1:2);
        
        gamma = sqrt(det(n12*C+G*N12));
        E12 = (G*M12-m12*S*C'*S')/gamma;
        F12 = (n12*C+G*N12)/gamma;
        
        CL=(M12*C+G*m12)*S*F12'*S';
        AL=E12*A*S*E12'*S';
        BL=F12*B*S*F12'*S';
    end

    function mask=setelems(mask,idx)
        mask(idx)=true;
    end

    function UP = BetatronPhaseUnwrap(P)
        % unwrap negative jumps in betatron
        %JUMPS = [0; diff(P)] < -1.e-5;
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
        catch %#ok<CTCH>
            phase=NaN(size(beta));
        end
        phase = BetatronPhaseUnwrap(phase);
    end

    function [beta,alpha,tune] = closure(AB)
        cosmu = (AB(1,1) + AB(2,2))/2;
        diff  = (AB(1,1) - AB(2,2))/2;
        sinmu = sign(AB(1,2))*sqrt(-AB(1,2)*AB(2,1)-diff*diff);
        alpha = (AB(1,1)-AB(2,2))/2/sinmu;
        beta = AB(1,2)/sinmu;
        tune = mod(atan2(sinmu,cosmu)/2/pi,1);
    end        

end
