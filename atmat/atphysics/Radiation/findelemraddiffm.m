function [B, M, r] = findelemraddiffm(ELEM,rin, varargin);
%FINDELEMRADDIFFM

MaxOrder = 2;

B = zeros(6);
M = eye(6);
r = rin;

switch ELEM.PassMethod
    case {'BndMPoleSymplectic4Pass', 'BndMPoleSymplectic4RadPass'}
        invrho = ELEM.BendingAngle/ELEM.Length;
        if isfield(ELEM,'T1')
            r = r + ELEM.T1';
        end
        if isfield(ELEM,'R1')
            r = ELEM.R1*r;
            B = ELEM.R1*B*ELEM.R1';
            M = ELEM.R1*M;
        end

        % Entrance edge effects
        E = eye(6);
        E(2,1) = invrho*tan(ELEM.EntranceAngle);
        if isfield(ELEM,'FullGap') & isfield(ELEM,'FringeInt1')
            E(4,3) = -invrho*tan(ELEM.EntranceAngle...
                -ELEM.FullGap*ELEM.FringeInt1...
                *(1+sin(ELEM.EntranceAngle)^2)/cos(ELEM.EntranceAngle)...
                /(1+r(5)));
        else
            E(4,3) = -invrho*tan(ELEM.EntranceAngle);
        end
        
        r = E*r;
        M = E*M;
        B = E*B*E';
        
        % Body of dipole
        [Bbody, Mbody, r] = findthickmpoleraddiffm(r,...
            ELEM.PolynomA, ELEM.PolynomB, ELEM.Length, invrho, ...
            ELEM.Energy, MaxOrder,ELEM.NumIntSteps);
        
        
        B = Mbody*B*Mbody'+Bbody;    
        M = Mbody*M;
        
        % Exit edge effects
        E = eye(6);
        E(2,1) = invrho*tan(ELEM.ExitAngle);
        if isfield(ELEM,'FullGap') & isfield(ELEM,'FringeInt2')
            E(4,3) = -invrho*tan(ELEM.ExitAngle...
                -ELEM.FullGap*ELEM.FringeInt1...
                *(1+sin(ELEM.ExitAngle)^2)/cos(ELEM.ExitAngle)...
                /(1+r(5)));
        else
            E(4,3) = -invrho*tan(ELEM.ExitAngle);
        end
        
        r = E*r;
        M = E*M;
        B = E*B*E';
        
        
        if isfield(ELEM,'R2')
            r = ELEM.R2*r;
            B = ELEM.R2*B*ELEM.R2';
            M = ELEM.R2*M;
        end
        
        if isfield(ELEM,'T2')
            r = r + ELEM.T2';
        end
        
    case {'StrMPoleSymplectic4Pass', 'StrMPoleSymplectic4RadPass'}

        if isfield(ELEM,'T1')
            r = r + ELEM.T1;
        end
        if isfield(ELEM,'R1')
            r = ELEM.R1*r;
            B = ELEM.R1*B*ELEM.R1';
            M = ELEM.R1*M;
        end

        % Body
        [Bbody, Mbody, r] = findthickmpoleraddifm(r,...
            ELEM.PolynomA, ELEM.PolynomB, ELEM.Length, invrho, ...
            ELEM.Energy, MaxOrder,ELEM.NumIntSteps);
        
        
        B = Mbody*B*Mbody'+Bbody;    
        M = Mbody*M;
        
        
        if isfield(ELEM,'R2')
            r = ELEM.R2*r;
            B = ELEM.R2*B*ELEM.R2';
            M = ELEM.R2*M;
        end
        
        if isfield(ELEM,'T2')
            r = r + ELEM.T2;
        end
        
    case 'BendLinearPass'
        % Add fields to element to make it 'BndMPoleSymplectic4Pass
        % compatible and call FINDELEMRADDIFFM recursively
        elem = ELEM;
        elem.PassMethod = 'BndMPoleSymplectic4Pass';
        elem.PolynomB = [0 0 0];
        elem.PolynomB = [0 elem.K 0];
        elem.NumIntSteps = 10;
        elem.MaxOrder = 2;
        [B, M, r] = findelemraddiffm(elem,rin);
    
    case 'IdentityPass'
        % Do nothing : [M, B, r] are the same at the exit 
    
    otherwise
        M = findelemm66(ELEM,ELEM.PassMethod,rin);
        r = feval(ELEM.PassMethod,ELEM,rin);
        B = M*B*M';
        
end