function varargout=atx(ring,varargin)
%ATX				computes and displays global information
%
%BEAMDATA=ATX(RING,DPP,REFPTS)
%
%RING:		AT structure
%DPP:		relative energy deviation (default: 0)
%REFPTS:    Index of elements (default: 1:length(ring))
%
%BEAMDATA is a MATLAB structure array with fields
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
%   M44         - 4x4 transfer matrix M from the beginning of RING
%                 to the entrance of the element for specified DP [2]
%   A           - 2x2 matrix A in [3]
%   B           - 2x2 matrix B in [3]
%   C           - 2x2 matrix C in [3]
%   gamma       - gamma parameter of the transformation to eigenmodes
%   mu          - [ mux, muy] horizontal and vertical betatron phase
%   beta        - [betax, betay] vector
%   alpha       - [alphax, alphay] vector
%
% From ohmienvelope:
%
%   beam66      - 6x6 equilibrium beam matrix
%   emit66      - 6x6 emittance projections on x and y + energy spread
%   beam44      - intersection of beam66 for dpp=0
%   emit44      - emittances of the projections of beam44 on x and y
%   modemit     - emittance of eigenmodes
%
%[BEAMDATA,PARAMS]=ATX(...)  Returns also a structure PM with fields
%   ll          - Circumference
%   alpha       - momentum compaction factor
%   nuh         - Tunes
%   nuv
%   fulltunes
%   fractunes
%   espread     - Energy spread
%   blength     - Bunch length
%   modemittance - Eigen emittances
%
%BEAMDATA=ATX(RING,DPP,REFPTS,RADRING,RADINDEX,CAVINDEX)
% Radiation must be turned on for emittance computation. This is done by
% default using the ATRADON function with default arguments. If this is not
% desired, this syntax allows to explicitly enter the radiative lattice
% with the indices of radiating elements and cavities
%
%
%BEAMDATA=ATX(RING,DPP,REFPTS,RADFUNCTION)
% RADFUNCTION is substituted to ATRADON to provide the radiative lattice
% and indices, in the form:
%        [RADRING,RADINDEX,CAVINDEX]=RADFUNCTION(RING)
%
% See also: ATLINOPT ATRADON OHMIENVELOPE

[energy,periods,voltage,~,eloss]=atenergy(ring);
[varargout{1:nargout}]=atx2(ring(:,1),energy,periods,voltage,eloss,varargin{:});
if nargout >= 1 && size(ring,2) > 1
    varargout{1}=repmat(varargout{1}(:),1,size(ring,2));
end

    function [linusr,pm]=atx2(ring,energy,periods,voltage,eloss,varargin)
        
        c=2.9987924e8;
        [dpp,refusr]=parseargs({0,1:length(ring)},varargin);
        if islogical(refusr)
            refusr(end+1,length(ring)+1)=false;
        else
            refusr=setelems(false(1,length(ring)+1),refusr);
        end
        refpts=setelems(refusr,[1 length(ring)+1]);
        keep=refusr(refpts);
        
        [lindata,tunes,xsi]=atlinopt(ring,dpp,refpts);
        coupled=max(abs(1-cat(1,lindata.gamma))) > 1.e-3; %#ok<*NOPRT>
        
        ts=periods*tunes;                       % fractional tunes
        fractunes=ts-fix(ts);
        
        circumference=periods*lindata(end).SPos;% circumference
        revperiod=lindata(end).SPos/c;
        
        transfer=lindata(end).M44;              % transfer matrix
        
        [beamA,beamB]=beam44(lindata(1));       % beam matrices
        
        closedorbit=lindata(1).ClosedOrbit';	% closed orbit
        
        dispersion=lindata(1).Dispersion';		% dispersion
        
        tuneper=lindata(end).mu/2/pi;			% tunes
        tunes=periods*tuneper;                  % tunes
        
        momcompact=mcf(ring);                   % momentum compaction
        
        chromaticity=xsi./tuneper;				% chromaticity
        
        synchrophase=asin(eloss/voltage);
        
        if nargout == 0
            display(coupled);
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
        
        if length(varargin)<3
            [ring2,radindex,cavindex]=atradon(ring);
        elseif isa(varargin{3},'function_handle')
            [ring2,radindex,cavindex]=varargin{3}(ring);
        else
            [ring2,radindex,cavindex]=deal(varargin{3:5});
        end
        
        if any(cavindex)
            try
                [envelope,espread,blength,m,T]=ohmienvelope(ring2,radindex,refpts);
                [chi,tns]=atdampingrates(m);
                fs=abs(tns(3))/revperiod;
                dampingtime=revperiod./chi;
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
                dampingtime=NaN(1,3);
                lindata=arrayfun(@deflt,lindata);
            end
        else
            blength=NaN;
            espread=NaN;
            fs=NaN;
            dampingtime=NaN(1,3);
            lindata=arrayfun(@deflt,lindata);
        end
        modemit=cat(1,lindata.modemit);
        modemittance=mean(modemit);
        modcoupling=mean(modemit(:,2)./modemit(:,1));
        projemit=cat(1,lindata.emit44);
        projemittance=projemit(1,:);
        projcoupling=mean(projemit(:,2)./projemit(:,1));
        if nargout==0
            display(energy);
            display(modemittance);
            display(modcoupling);
            display(projemittance);
            display(projcoupling);
            display(dampingtime);
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
                'dampingtime',dampingtime,...
                'espread',espread,...
                'blength',blength,...
                'modemittance',modemittance,...
                'energy',energy,...
                'fs',fs,...
                'eloss',eloss,...
                'synchrophase',synchrophase);
        end
        
        function lind=process(bm66,T,lind)
            if coupled             % bm44 in the intersection of beam66
                siginv=inv(bm66);   % with dp/p==0  ( 4x4 betatron emittance)
                bm44=inv(siginv(1:4,1:4));
            else
                siginv=inv(bm66([1 2 5 6],[1 2 5 6]));
                bm44=[inv(siginv(1:2,1:2)) zeros(2,2);zeros(2,4)];
            end
            % [beamA,beamB]=beam44(lindata(i));  % betatron 4x4 coupled matrices
            %   beamA       - 2x2 beam matrix for mode A
            %   beamB       - 2x2 beam matrix for mode B
            % Eigen emittances: Solve bm44 = modemit(1)*beamA + modemit(2)*beamB;
            % lindata(i).modemit=([beamA(:) beamB(:)]\bm44(:))';
            % aa=amat(T(:,:,i)*m*inv(T(:,:,i))); %#ok<MINV>
            aa=amat(T*m*jmt'*T'*jmt);
            nn=-aa'*jmt*bm66*jmt*aa;
            lind.modemit=0.5*[nn(1,1)+nn(2,2) nn(3,3)+nn(4,4) nn(5,5)+nn(6,6)];
            lind.beam66=bm66;
            lind.beam44=bm44;
            % Projected betatron emittances
            lind.emit44=sqrt([det(bm44(1:2,1:2)) det(bm44(3:4,3:4))]);
            % Projected full emittances (energy spread included)
            lind.emit66=sqrt([det(bm66(1:2,1:2)) det(bm66(3:4,3:4)) det(bm66(5:6,5:6))]);
        end
        
        function lind=deflt(lind)
            lind.modemit=NaN(1,3);
            lind.emit44=NaN(1,2);
        end
        
        function [chi,nu]=atdampingrates(m66)
            %find tunes and damping rates from one map matrix with radiation
            aa=amat(m66);
            
            Rmat=aa\m66*aa;
            
            [chi,nu]=cellfun(@decode,num2cell(reshape(1:size(m66,1),2,[]),1));
            
            function [chi,nu]=decode(range)
                matr=Rmat(range,range);
                nu=atan2(matr(1,2)-matr(2,1),matr(1,1)+matr(2,2))/2/pi;
                chi=-log(sqrt(det(matr)));
            end
        end
        
        function mask=setelems(mask,idx)
            mask(idx)=true;
        end
    end
end


