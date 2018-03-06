function rerr = atseterrortable(r,ErrorTable,varargin)
% ATSETERRORTABLE Sets errors listed in table ErrorsTable in the AT lattice r
%
% rerr = ATSETERRORTABLE(r,ErrorTable,varargin)
%
% INPUTS:
% r : AT lattice (cell array of structures)
% ErrorTable: matlab table with the columns described below:
%
%    FamNames    Index    X    Y    S    Roll    Pitch    Yaw    DAngle_Angle    Bend_Roll        DK_K         b_n_systematic    a_n_systematic     b_n_random       a_n_random      BPM_Offset    BPM_Gain    BPM_Reading    BPM_Rotation
%    ________    _____    _    _    _    ____    _____    ___    ____________    _________    _____________    ______________    ______________    _____________    _____________    __________    ________    ___________    ____________
%    SDHI         1       0    0    0    0       0        0      0               0            [1x20 double]    [1x20 double]     [1x20 double]     [1x20 double]    [1x20 double]    0    0        0    0      0    0         0
%    BPM          2       0    0    0    0       0        0      0               0            [1x20 double]    [1x20 double]     [1x20 double]     [1x20 double]    [1x20 double]    0    0        0    0      0    0         0
%    SD1U         3       0    0    0    0       0        0      0               0            [1x20 double]    [1x20 double]     [1x20 double]     [1x20 double]    [1x20 double]    0    0        0    0      0    0         0
%    S4           4       0    0    0    0       0        0      0               0            [1x20 double]    [1x20 double]     [1x20 double]     [1x20 double]    [1x20 double]    0    0        0    0      0    0         0
%    SD2U         5       0    0    0    0       0        0      0               0            [1x20 double]    [1x20 double]     [1x20 double]     [1x20 double]    [1x20 double]    0    0        0    0      0    0         0
%    QF2          6       0    0    0    0       0        0      0               0            [1x20 double]    [1x20 double]     [1x20 double]     [1x20 double]    [1x20 double]    0    0        0    0      0    0         0
%    SD2D         7       0    0    0    0       0        0      0               0            [1x20 double]    [1x20 double]     [1x20 double]     [1x20 double]    [1x20 double]    0    0        0    0      0    0         0
%    BPM          8       0    0    0    0       0        0      0               0            [1x20 double]    [1x20 double]     [1x20 double]     [1x20 double]    [1x20 double]    0    0        0    0      0    0         0
%    SDA          9       0    0    0    0       0        0      0               0            [1x20 double]    [1x20 double]     [1x20 double]     [1x20 double]    [1x20 double]    0    0        0    0      0    0         0
%    S6          10       0    0    0    0       0        0      0               0            [1x20 double]    [1x20 double]     [1x20 double]     [1x20 double]    [1x20 double]    0    0        0    0      0    0         0
%
% where:
% FamName: reported for simpler reading of the table
% Index: index in the r AT lattice, used to assign the errors
% X     is the radial displacement in [m] (T1, T2 fields are modified)
% Y     is the vertical displacement in [m] (T1, T2 fields are modified)
% S     is the longitudinal displacement in [m] (drift lengths and T1 T2 are modified)
% Roll  is the rotation about s in [m] (R1,R2 fields are modified)
% Pitch is the rotation about x in [m] (T1, T2, R1, R2 fields are modified)
%       Pitch depends on length of magnets. 
%       This function does not take care of sliced magnets.        
% Yaw   is the rotation about y in [m] (T1, T2, R1, R2 fields are modified)
%       Yaw depends on length of magnets. 
%       This function does not take care of sliced magnets.        
% DAngle_Angle is the BendingAngle error in [rad] respect to BendingAngle
%       if existing (assigned as Delta PolynomB(1))
%       Angle depends on length of magnets. 
%       This function does not take care of sliced magnets.        
% Bend_Roll is the Bending magnet rotation [rad] implemented using
%       polynomB, polynomA (no change of reference trajectory, lattice survey unchanged)
% DK_K  is the field integral error in [1/m,1/m2,1/m3...] respect to PolynomB
% b_n_systematic are normal systematic multipole errors summed to PolynomB
% a_n_systematic are skew systematic multipole errors summed to PolynomA
% b_n_random     are normal random multipole errors summed to PolynomB
% a_n_random     are skew random multipole errors summed to PolynomA
% BPM_Offset    are BPM offsets
% BPM_Gain      are BPM scale errors
% BPM_Reading   are BPM random reading errors (sigma)
% BPM_Rotation  are BPM rotations
%
%OPTIONAL INPUT:
% verbose: print what is happening
%
% NOTE: Errors are set as in the table. if slice magnet elements or magnets on girders need to
%       have the same errors this should already be taken in account
%       in the input table.
%
% NOTE: Pitch and Yaw are implementd using T1 =! T2
%
% see also: atsetshift atsetbpmerr atset_s_shift atsettilt atsetpitch
%          atsetyaw

p = inputParser;
addRequired(p,'r',@iscell);
addOptional(p,'verbose',false,@islogical);

parse(p,r,varargin{:});

r=p.Results.r;
verbose=p.Results.verbose;


rerr = r;        % initialize lattice with errors

% set BPM errors
if find([ErrorTable.BPM_Offset ErrorTable.BPM_Reading ErrorTable.BPM_Rotation])
    rerr = atsetbpmerr(rerr,ErrorTable.Index,...
        ErrorTable.BPM_Offset(:,1),...
        ErrorTable.BPM_Offset(:,2),...
        ErrorTable.BPM_Gain(:,1),...
        ErrorTable.BPM_Gain(:,2),...
        ErrorTable.BPM_Reading(:,1),...
        ErrorTable.BPM_Reading(:,2),...
        ErrorTable.BPM_Rotation);
    if verbose, disp('set BPM errors'); end
end

% set X,Y shifts
if find([ErrorTable.X ErrorTable.Y])
    rerr = atsetshift(rerr,ErrorTable.Index,ErrorTable.X,ErrorTable.Y);
    if verbose, disp('set Dx Dy alignment errors'); end
end

% set s displacement
if find(ErrorTable.S)
    rerr = atset_s_shift(rerr,ErrorTable.Index,ErrorTable.S);
    if verbose, disp('set Ds alignment errors'); end
end

% set S rotations
if find(ErrorTable.Roll)
    rerr = atsettilt(rerr,ErrorTable.Index,ErrorTable.Roll);
    if verbose, disp('set roll errors (rotation about s)'); end
end

% set pitch
if find(ErrorTable.Pitch)
    rerr = atsetpitch(rerr,ErrorTable.Index,ErrorTable.Pitch);
    if verbose, disp('set pitch errors (rotation about x (radial))'); end
end

% set yaw
if find(ErrorTable.Yaw)
    rerr = atsetyaw(rerr,ErrorTable.Index,ErrorTable.Yaw);
    if verbose, disp('set Yaw errors (rotation about y (vertical))'); end
end

% set yaw
if find(ErrorTable.Yaw)
    for im = 1:length(ErrorTable.Index)
        if isfield(r{im},'BendingAngle')
            BA = r{im}.BendingAngle;
            L = r{im}.Length;
            ErrorTable.b_n_systematic(im,1) = ...
                ErrorTable.b_n_systematic(im,1) + ...
                ErrorTable.DAngle_Angle(im,1) * BA / L;
        end
    end
    
    if verbose, disp('converted DAngle_Angle to b_n_syst(:,1)'); end
end

% set gradient errors
if find(ErrorTable.DK_K)
    
    magind = ErrorTable.Index;
    
    for im = 1:length(magind)
        if isfield(r{magind(im)},'PolynomB') % if magnet
            % get DK_K errors
            DK_K   = ErrorTable.DK_K(magind(im),:);
            % find last non zero
            lastnonzero = find(DK_K,1,'last');
            if ~isempty(lastnonzero) % if there are gradient errors
                % set maxorder
                magmaxord = rerr{magind(im)}.MaxOrder;
                maxord = max(lastnonzero,magmaxord) -1;
                rerr=atsetfieldvalues(rerr,magind(im),'MaxOrder',maxord); % get max order
                % pad polynomA polynomB with zeros
                rerr(magind(im))=PadPolynomAB(rerr(magind(im)),maxord);
                % get existing polynomB
                pb0 = atgetfieldvalues(r,magind(im),'PolynomB',{1,1:maxord+1});
                if isnumeric(pb0),pb0={pb0}; end;
                % multiply by fractional errors
                newkl = pb0{1}.*(1+DK_K(1:(maxord+1)));
                % set new PolynomB with errors
                rerr = atsetfieldvalues(rerr,magind(im),'PolynomB',{1,1:maxord+1},newkl);
                % disp(pb0{1}-newkl)
            end
        end
    end
    if verbose, disp('set gradient errors'); end
end

% set systematic multipole errors
if find([ErrorTable.b_n_systematic ErrorTable.b_n_random])
    
    magind = ErrorTable.Index;
    
    for im = 1:length(magind)
        if isfield(r{magind(im)},'PolynomB') % if magnet
            % get normal multipoles errors
            b_n_syst   = ErrorTable.b_n_systematic(magind(im),:);
            a_n_syst   = ErrorTable.a_n_systematic(magind(im),:);
            % find last non zero
            lb = find(b_n_syst,1,'last'); if isempty(lb), lb = 0; end
            la = find(a_n_syst,1,'last'); if isempty(la), la = 0; end
            lastnonzero = max(la,lb);
            if ~isempty(lastnonzero) % if there are multipole errors
                % set maxorder
                maxord = lastnonzero -1;
                rerr=atsetfieldvalues(rerr,magind(im),'MaxOrder',maxord); % get max order
                % pad polynomA polynomB with zeros
                rerr(magind(im)) = PadPolynomAB(rerr(magind(im)),maxord+1);
                % get existing polynomB
                pb0 = atgetfieldvalues(rerr,magind(im),'PolynomB',{1,1:maxord+1});
                pa0 = atgetfieldvalues(rerr,magind(im),'PolynomA',{1,1:maxord+1});
                % add errors
                newb = pb0{1} + b_n_syst(1:(maxord+1));
                newa = pa0{1} + a_n_syst(1:(maxord+1));
                
                % set new PolynomB with errors
                rerr = atsetfieldvalues(rerr,magind(im),'PolynomB',{1,1:maxord+1},newb);
                rerr = atsetfieldvalues(rerr,magind(im),'PolynomA',{1,1:maxord+1},newa);
                
            end
        end
    end
    if verbose, disp('set systematic multipole errors'); end
end


% set random multipole errors
if find([ErrorTable.b_n_random ErrorTable.a_n_random])
    
    magind = ErrorTable.Index;
    
    for im = 1:length(magind)
        if isfield(r{magind(im)},'PolynomB') % if magnet
            % get normal multipoles errors
            b_n_rand   = ErrorTable.b_n_random(magind(im),:);
            a_n_rand   = ErrorTable.a_n_random(magind(im),:);
            % find last non zero
            lb = find(b_n_syst,1,'last'); if isempty(lb), lb = 0; end
            la = find(a_n_syst,1,'last'); if isempty(la), la = 0; end
            lastnonzero = max(la,lb);
            if ~isempty(lastnonzero) % if there are multipole errors
                % set maxorder
                maxord = lastnonzero -1;
                rerr=atsetfieldvalues(rerr,magind(im),'MaxOrder',maxord); % get max order
                % pad polynomA polynomB with zeros
                rerr(magind(im)) = PadPolynomAB(rerr(magind(im)),maxord+1);
                % get existing polynomB
                pb0 = atgetfieldvalues(rerr,magind(im),'PolynomB',{1,1:maxord+1});
                pa0 = atgetfieldvalues(rerr,magind(im),'PolynomA',{1,1:maxord+1});
                % add errors
                newb = pb0{1} + b_n_rand(1:(maxord+1));
                newa = pa0{1} + a_n_rand(1:(maxord+1));
                
                % set new PolynomB with errors
                rerr = atsetfieldvalues(rerr,magind(im),'PolynomB',{1,1:maxord+1},newb);
                rerr = atsetfieldvalues(rerr,magind(im),'PolynomA',{1,1:maxord+1},newa);
                % disp(pb0{1}-newkl)
            end
        end
    end
    if verbose, disp('set random multipole errors'); end
end

% set Dipole S rotations (modify polynomB polynomA, must rotate all multipoles  also.)
if find(ErrorTable.Bend_Roll)
    rerr = atsettiltdipole(rerr,ErrorTable.Index,ErrorTable.Bend_Roll);
    if verbose, disp('set Bending magnets roll errors (rotation about s) without changing reference trajectory'); end
end

end



function r=PadPolynomAB(r,maxord)
% pads PolynomB and A to have the same number of elements


r=cellfun(@(a)padpol(a,maxord),r,'un',0);

end

function a=padpol(a,maxord)

if isfield(a,'PolynomB')
    try
        lpa=length(a.PolynomA);
        lpb=length(a.PolynomB);
    catch
        a.PolynomA=[0];
        a.PolynomB=[0];
        a.MaxOrder=0;
        a.NumIntSteps=1;
        lpa=length(a.PolynomA);
        lpb=length(a.PolynomB);
    end
    
    if lpa<lpb
        a.PolynomA=[a.PolynomA,zeros(1,lpb-lpa)];
    elseif lpa>lpb
        a.PolynomB=[a.PolynomB,zeros(1,lpa-lpb)];
    end
    
    padlen = length(a.PolynomB);
    if maxord > padlen
        a.PolynomB=[a.PolynomB,zeros(1,maxord-padlen)];
        a.PolynomA=[a.PolynomA,zeros(1,maxord-padlen)];
    end
    
end

end

