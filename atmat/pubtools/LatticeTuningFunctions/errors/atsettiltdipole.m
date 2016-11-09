function ring=atsettiltdipole(varargin)
%ATSETTILTDIPOLE sets the entrance and exit rotation matrices
% of an element or a group of elements in THERING
%
% RING=ATSETTILTDIPOLE(RING,ELEMINDEX, PSI)
% ELEMINDEX contains indexes of elements to be rotated
% PSI - angle(s) of rotation in RADIANS
%   POSITIVE PSI corresponds to a CORKSCREW (right)
%   rotation of the ELEMENT looking in the direction of the beam.
%   (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
%   The misalgnment matrixes are stored in fields R1 and R2
%   R1 = [  cos(PSI) sin(PSI); -sin(PSI) cos(PSI) ]
%   R2 = R1'
%
% the rotated dipole gives a kick in the vertical and horizontal plane.
%
% ATSETTILTDIPOLE(ELEMINDEX, PSI) Uses the global variable THERING
%
% See also ATSETSHIFT ATSETTILT

global THERING
if ~iscell(varargin{1})
    THERING=atsettilt(THERING,varargin{:});
else
    [ring,idx,rot]=deal(varargin{:});
    
    if length(rot) == 1
        rot=rot*ones(size(idx));
    elseif length(rot) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(rot))
    end
    
    %     tic;
    %     for i = 1:length(idx)
    %         ring{idx(i)}=attiltelemdip(ring{idx(i)},rot(i));
    %     end
    %     toc
    %tic;
    ring(idx)=cellfun(@(el,rot)attiltelemdip(el,rot),ring(idx),num2cell(rot(:)),'un',0);
    %toc
    
end

end

function elem = attiltelemdip(elem,rots)
%ATTILTELEMdip set new rotation parameters
%NEWELEM=ATTILTELEMdip(OLDELEM,ROTS)
%
% ROTS - rotation angle in RADIANS
%   POSITIVE ROTS corresponds to a CORKSCREW (right)
%   rotation of the ELEMENT looking in the direction of the beam.
%   (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
%   The rotation matrixes are stored in fields R1 and R2
%   R1 = [  cos(PSI) sin(PSI); -sin(PSI) cos(PSI) ]
%   R2 = R1'
%See also: atshiftelem, atmodelem

if ~isfield(elem,'BendingAngle')% rotate reference
     
    C=cos(rots);
    S=sin(rots);
    RM = diag([C C C C 1 1]);
    RM(1,3) = S;
    RM(2,4) = S;
    RM(3,1) = -S;
    RM(4,2) = -S;
    elem.R1=RM;
    elem.R2=RM';
    
else% rotate multipoles
    
    % bending angle
    bb=-elem.BendingAngle/elem.Length;
    
    % horizontal kick
    elem.('PolynomB')(1)=-(1-cos(rots)).*bb;
    % vertical kick
    elem.('PolynomA')(1)=sin(rots).*bb;
    
    % rotate all other multipole components ( combined function magnets)
    Lpb=length(elem.('PolynomB'));
    
    elem=padpol(elem);
    
    rotm=-rots*[2:Lpb];
    
    elem.('PolynomB')(2:end)=cos(rotm).*elem.('PolynomB')(2:end)-sin(rotm).*elem.('PolynomA')(2:end);
    elem.('PolynomA')(2:end)=sin(rotm).*elem.('PolynomB')(2:end)+cos(rotm).*elem.('PolynomA')(2:end);
    
end

elem.RotAboutS=rots;

end



function a=padpol(a)

if isfield(a,'PolynomB')
    
    lpa=length(a.PolynomA);
    lpb=length(a.PolynomB);
    
    if lpa<lpb
        a.PolynomA=[a.PolynomA,zeros(1,lpb-lpa)];
    elseif lpa>lpb
        a.PolynomB=[a.PolynomB,zeros(1,lpa-lpb)];
    end
    
end
end
