function [lindt,pm]=atx(ring,varargin)
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
%   modemit     - [emitA emitB] emittance of modes A and B (should be constant)
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
% See also: ATREADBETA ATLINOPT OHMIENVELOPE ATMODUL

if isvector(ring)
    [lindt,pm]=atx2(ring,nargout,1,varargin{:});
else
    periods=size(ring,2);
    ring1=ring(:,1);
    cavindex=findcells(ring1,'HarmNumber');
    ring2=setcellstruct(ring1,'HarmNumber',cavindex,getcellstruct(ring1,'HarmNumber',cavindex)/periods);
    [lindt,pm]=atx2(ring2,nargout,periods,varargin{:});
    lindt=repmat(lindt(:),1,periods);
end

function [linusr,pm]=atx2(ring,outargs,periods,dpp,refusr)

if nargin < 5
    refpts=1:length(ring)+1;
    keep=1:length(ring);
else
    [refpts,dum,keep]=unique([1 refusr length(ring)+1]);
    keep=keep(2:end-1);
end
if nargin < 4, dpp=0; end

[lindata,tunes,xsi]=atlinopt(ring,dpp,refpts);
%coupled=std(cat(1,lindata.gamma),1) > 1.e-5
coupled=max(abs(1-cat(1,lindata.gamma))) > 1.e-5; %#ok<*NOPRT>

ts=periods*tunes;						% fractional tunes
fractunes=ts-fix(ts);

circumference=periods*lindata(end).SPos;			% circumference

transfer=lindata(end).M44;			% transfer matrix

[beamA,beamB]=beam44(lindata(1));

closedorbit=lindata(1).ClosedOrbit';		% closed orbit

dispersion=lindata(1).Dispersion';		% dispersion

tuneper=lindata(end).mu/2/pi;			% tunes
tunes=periods*tuneper;			% tunes

momcompact=mcf(ring);				% momentum compaction

chromaticity=xsi./tuneper;				% chromaticity

if outargs == 0
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
end

[ring2,radindex,cavindex]=atradon(ring);

if ~isempty(cavindex)
    try
        [envelope,espread,blength]=ohmienvelope(ring2,radindex,refpts);
        for i=1:length(lindata)
            [beamA,beamB]=beam44(lindata(i));  % betatron 4x4 coupled matrices
            %   beamA       - 2x2 beam matrix for mode A
            %   beamB       - 2x2 beam matrix for mode B
            bm66=envelope(i).R;    % bm66 is the 6x6 beam matrix
            %         coupled=true;
            if coupled             % bm44 in the intersection of beam66
                siginv=inv(bm66);   % with dp/p==0  ( 4x4 betatron emittance)
                bm44=inv(siginv(1:4,1:4));
            else
                siginv=inv(bm66([1 2 5 6],[1 2 5 6]));
                bm44=[inv(siginv(1:2,1:2)) zeros(2,2);zeros(2,4)];
            end
            lindata(i).beam66=bm66;
            lindata(i).beam44=bm44;
            % Eigen emittances: Solve bm44 = modemit(1)*beamA + modemit(2)*beamB;
            lindata(i).modemit=([beamA(:) beamB(:)]\bm44(:))';
            % Projected betatron emittances
            lindata(i).emit44=sqrt([det(bm44(1:2,1:2)) det(bm44(3:4,3:4))]);
            % Projected full emittances (energy spread included)
            lindata(i).emit66=sqrt([det(bm66(1:2,1:2)) det(bm66(3:4,3:4)) bm66(5,5)]);
            lindata(i).tilt=envelope(i).Tilt;
        end
        modemit=cat(1,lindata.modemit);
        %emit6=cat(1,lindata.emit66);
        modemittance=mean(modemit);
        modcoupling=mean(modemit(:,2)./modemit(:,1));
        projemit=cat(1,lindata.emit44);
        projemittance=projemit(1,:);
        projcoupling=mean(projemit(:,2)./projemit(:,1));
        if outargs==0
            display(modemittance);
            display(modcoupling);
            display(projemittance);
            display(projcoupling);
        end
    catch err
        warning('atx:unstable','Emittance computation failed:\n%s',err.message);
        blength=NaN;
        espread=NaN;
        modemittance=NaN(1,2);
        for i=1:length(lindata)
            %		 lindata(i).beam66=bm66;
            %		 lindata(i).beam44=bm44;
            lindata(i).modemit=NaN(2,1);
            lindata(i).emit44=NaN(2,1);
            %		 lindata(i).emit66=sqrt([det(bm66(1:2,1:2)) det(bm66(3:4,3:4)) bm66(5,5)]);
            %        lindata(i).tilt=envelope(i).Tilt;
        end
    end
else
    blength=NaN;
    espread=NaN;
    modemittance=NaN(1,2);
    for i=1:length(lindata)
        %		 lindata(i).beam66=bm66;
        %		 lindata(i).beam44=bm44;
        lindata(i).modemit=NaN(2,1);
        lindata(i).emit44=NaN(2,1);
        %		 lindata(i).emit66=sqrt([det(bm66(1:2,1:2)) det(bm66(3:4,3:4)) bm66(5,5)]);
        %        lindata(i).tilt=envelope(i).Tilt;
    end
end
linusr=lindata(keep);
pm=struct('ll',circumference,'alpha',momcompact,...
    'fractunes',fractunes,...
    'fulltunes',tunes,...
    'nuh',tunes(1),'nuv',tunes(2),...
    'espread',espread,...
    'blength',blength,...
    'modemittance',modemittance);
if outargs==0
    display(espread);
    display(blength);
end


