function varargout=atx(ring,varargin)
%ATX Computes and displays global information
%
%BEAMDATA=ATX(RING) Computes linear optics, equilibrium emittances and damping times
%
% RING:     Ring description. If RING is 6D, it is used in OHMIENVELOPE and
%           a 4D copy may be used for linear optics computation.
%           If RING is 4D, a 6D copy with default options is used for
%           OHMIENVELOPE
%
%BEAMDATA=ATX(RING,DP,REFPTS)
%BEAMDATA=ATX(RING,REFPTS)
%   Specify the points of interest (Default: 1:length(RING)+1)
%
%BEAMDATA=ATX(RING,DP,...) 
%BEAMDATA=ATX(RING,...,'dp',DPP)
%   Specify the momentum deviation (Default: 0)
%
%BEAMDATA=ATX(RING,...,'dct',DCT)
%   Specify the path lengthening
%
%BEAMDATA=ATX(RING,...,'method',OPTICSFUN)
%   Specify the method for linear optics. Default: @atlinopt6
%   Allowed values are @atlinopt2, @atlinopt4, @atlinopt6
%
%BEAMDATA=ATX(RING,...,'6d')
%   By default, linear optics is computed with the 4d version of the lattice.
%   If method is @atlinopt6 (the default), when specifying '6d' the optics
%   is computed from the 6d version of the lattice. 
%
%ELEMDATA is a MATLAB structure array as long as the numver of refpoints
%with fields:
%
% From atlinopt:
%
%   ElemIndex   - ordinal position in the RING
%   SPos        - longitudinal position [m]
%   ClosedOrbit - closed orbit column vector with
%                 components x, px, y, py (momentums, NOT angles)
%   Dispersion  - dispersion orbit position vector with
%                 components eta_x, eta_prime_x, eta_y, eta_prime_y
%                 calculated with respect to the closed orbit with
%                 momentum deviation DP
%   M           - 4x4 transfer matrix M from the beginning of RING
%                 to the entrance of the element for specified DP [2]
%   mu          - [ mux, muy] horizontal and vertical betatron phase
%   beta        - [betax, betay] vector
%   alpha       - [alphax, alphay] vector
%
%   Other fields depend on the selected optics method, @atlinopt4 or
%   @atlinopt6
%
% From ohmienvelope:
%
%   beam66      - 6x6 equilibrium beam matrix
%   emit66      - 6x6 emittance projections on x-x', y-y' and energy spread
%   beam44      - intersection of beam66 with dpp=0 (monochromatic beam)
%   emit44      - emittances of the projections of beam44 on x and y
%   modemit     - emittance of the 3 eigenmodes
%
%[ELEMDATA,RINGDATA]=ATX(...)  Returns also a structure RINGDATA
%with fields:
%
%   ll          - Circumference
%   alpha       - momentum compaction factor
%   fractunes
%   fulltunes
%   nuh         - Tunes
%   nuv
%   chromaticity
%   dampingtime
%   espread     - Energy spread
%   blength     - Bunch length
%   energy
%   fs          - synchrotron frequency
%   eloss       - energy loss/turn
%   synchrophase- synchronous phase
%   modemittance- Eigen emittances
%   momcompact  - momentum compaction factor
%
%
%The following options are kept for backwards compatibility but are
%deprecated:
%
%BEAMDATA=ATX(RING,DPP,REFPTS,RADRING,RADINDEX,CAVINDEX)
% Radiation must be turned on for emittance computation. This is done by
% default using the ATRADON function with default arguments. If this is not
% desired, this syntax allows to explicitly enter the radiative lattice
%
%BEAMDATA=ATX(RING,DPP,REFPTS,RADFUNCTION)
% RADFUNCTION is substituted to ATRADON to provide the radiative lattice
% and indices, in the form:
%        [RADRING,RADINDEX,CAVINDEX]=RADFUNCTION(RING)
%
% See also atlinopt atradon ohmienvelope ringpara atsummary

[go6d,vargs]=getflag(varargin,'6d');
[linopt,vargs]=getoption(vargs,'method',@atlinopt6);
vargs=getdparg(vargs);
[refusr,vargs]=getargs(vargs,1:length(ring),'check',@(arg) isnumeric(arg) || islogical(arg));
if ~any(cellfun(@(m) isequal(m,linopt), {@atlinopt2,@atlinopt4,@atlinopt6}))
    error('AT:WrongArgument','''method'' must be @atlinopt2, @atlinopt4 or @atlinopt6')
end
if go6d && ~isequal(linopt,@atlinopt6)
    error('AT:WrongArgument','''6d'' needs the @atlinopt6 method')
end
if ~isempty(vargs) && isa(vargs{1},'function_handle')
    [varargout{1:nargout}] = go(ring,vargs{1}(ring),refusr,vargs{2:end});
elseif ~isempty(vargs) && iscell(vargs{1})
    [varargout{1:nargout}] = go(ring,vargs{1},refusr,vargs{4:end});
else
    if atGetRingProperties(ring,'is_6d')
        [varargout{1:nargout}] = go(atdisable_6d(ring),ring,refusr,vargs{:});
    else
        [varargout{1:nargout}] = go(ring,atenable_6d(ring),refusr,vargs{:});
    end
end

    function [linusr,pm]=go(roff,ron,refusr,varargin)
        [dpargs,varargs]=getoption(varargin,{'dp','dct','df'});
        [cavargs,varargs]=getoption(varargs,{'cavpts'});
        if ~isempty(dpargs)
            ron=atsetcavity(ron,'Frequency','nominal',dpargs{:},cavargs{:});
        end
        [energy,periods,has_cavity,voltage,cell_length,cell_revfreq]=...
            atGetRingProperties(ron,'Energy','Periodicity','has_cavity',...
            'rf_voltage','cell_length','cell_revolution_frequency');
        radindex=atgetcells(ron,'PassMethod',@(elem,pass) endsWith(pass,'RadPass'));
        eloss=atgetU0(ron,'method','tracking');
        % Add 1 and length(ring)+1 to refpts
        if islogical(refusr)
            refusr(end+1,length(roff)+1)=false;
        else
            refusr=setelems(false(1,length(roff)+1),refusr);
        end
        refpts=setelems(refusr,[1 length(roff)+1]);
        keep=refusr(refpts);
        
        orbit4 = findorbit(roff,dpargs{:});
        dpp=orbit4(5);
        if go6d
            [ringdata,lindata]=linopt(ron,refpts,varargs{:},'get_chrom');
        else
            [ringdata,lindata]=linopt(roff,refpts,'orbit',orbit4,dpargs{:},varargs{:},'get_chrom');
        end
        
        ts=periods*ringdata.tune;               % fractional tunes
        fractunes=ts-fix(ts);
        
        circumference=periods*cell_length;      % circumference
        
        transfer=lindata(end).M;                % transfer matrix
        
        if isequal(linopt,@atlinopt6)           % beam matrices
            beamA=lindata(1).R(:,:,1);
            beamB=lindata(1).R(:,:,2);
        else
            [beamA,beamB]=beam44(lindata(1));
        end
        
        closedorbit=lindata(1).ClosedOrbit';	% closed orbit
        
        dispersion=lindata(1).Dispersion';		% dispersion
        
        tuneper=lindata(end).mu/2/pi;
        tunes=periods*tuneper;                  % tunes with integer part
        
        momcompact=mcf(roff,dpp);               % momentum compaction
        
        chromaticity=[periods*ringdata.chromaticity;ringdata.chromaticity./tuneper];		% chromaticity
        
        synchrophase=asin(eloss/voltage);       % synchronous phase
        
        if nargout == 0
            display(fractunes);
            display(circumference);
            display(transfer);
            display(beamA);
            display(beamB);
            display(closedorbit);
            display(dispersion);
            display(tuneper);
            display(tunes);
            display(momcompact);
            display(chromaticity);
            display(eloss);
            display(synchrophase);
        end
        
        if has_cavity
            try
                [envelope,espread,blength,m,T]=ohmienvelope(ron,radindex,refpts);
                [tns,chi]=atdampingrates(m);
                fs=abs(tns(3))*cell_revfreq;
                alpha=chi*cell_revfreq;
                dampingtime=1.0./alpha;
                dampingJ=4.0*alpha/sum(alpha);
                if any(radindex)
                    jmt=jmat(3);
                    lindata=cellfun(@process,{envelope.R},reshape(num2cell(T,[1 2]),1,[]),num2cell(lindata));
                else
                    lindata=arrayfun(@deflt,lindata);
                end
            catch err
                warning('atx:unstable','Emittance computation failed:\n%s\n',...
                    err.message);
                blength=NaN;
                espread=NaN;
                fs=NaN;
                dampingtime=NaN(1,3);
                dampingJ=NaN(1,3);
                lindata=arrayfun(@deflt,lindata);
            end
        else
            blength=NaN;
            espread=NaN;
            fs=NaN;
            dampingtime=NaN(1,3);
            dampingJ=NaN(1,3);
            lindata=arrayfun(@deflt,lindata);
        end
        % Average the mode emittance along the ring
        modemit=cat(1,lindata.modemit);
        modemittance=mean(modemit);
        modcoupling=modemittance(2)./modemittance(1);
        % Emittance projections on x-x', z-z' and delta-ctau
        projemit=cat(1,lindata.emit66);
        projemittance=projemit(1,:);
        % Averaged X-Z coupling (fluctuating...)
        projcoupling=mean(projemit(:,2)./projemit(:,1));
        if nargout==0
            display(energy);
            display(modemittance);
            display(modcoupling);
            display(projemittance);
            display(projcoupling);
            display(dampingtime);
            display(dampingJ);
            display(fs);
            display(espread);
            display(blength);
        end
        if nargout>=1
            linusr=lindata(keep);
        end
        if nargout>=2
            pm=struct('ll',circumference,'alpha',momcompact,...
                'fractunes',fractunes,...
                'fulltunes',tunes,...
                'nuh',tunes(1),'nuv',tunes(2),...
                'chromaticity',chromaticity(1,:),...
                'dampingtime',dampingtime,...
                'dampingJ',dampingJ,...
                'espread',espread,...
                'blength',blength,...
                'modemittance',modemittance,...
                'energy',energy,...
                'fs',fs,...
                'eloss',eloss,...
                'synchrophase',synchrophase,...
                'momcompact',momcompact);
        end
        
        function lind=process(bm66,T,lind)            
            % aa=amat(T*m*jmt'*T'*jmt);     % symplectic only
            aa=amat(T*m/T);                 % even for non-symplectic
            nn=-aa'*jmt*bm66*jmt*aa;
            % Mode emittances: should be constant
            memit=0.5*[nn(1,1)+nn(2,2) nn(3,3)+nn(4,4) nn(5,5)+nn(6,6)];
            lind.modemit=memit;
            if memit(2)/memit(1) > 1.e-4
                siginv=inv(bm66);   % with dp/p==0  ( 4x4 betatron emittance)
                bm44=inv(siginv(1:4,1:4));
            else
                siginv=inv(bm66([1 2 5 6],[1 2 5 6]));
                bm44=[inv(siginv(1:2,1:2)) zeros(2,2);zeros(2,4)];
            end
            % 6x6 full beam matrix
            lind.beam66=bm66;
            % 4x4 monochromatic beam matrix
            lind.beam44=bm44;
            % Monochromatic betatron emittances
            lind.emit44=sqrt([det(bm44(1:2,1:2)) det(bm44(3:4,3:4))]);
            % Projected full emittances (energy spread included)
            lind.emit66=sqrt([det(bm66(1:2,1:2)) det(bm66(3:4,3:4)) det(bm66(5:6,5:6))]);
        end
        
        function lind=deflt(lind)
            lind.modemit=NaN(1,3);
            lind.emit44=NaN(1,2);
            lind.emit66=NaN(1,3);
        end
        
        function mask=setelems(mask,idx)
            mask(idx)=true;
        end
    end
end
