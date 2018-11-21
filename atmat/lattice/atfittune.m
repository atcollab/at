function newring=atfittune(ring,varargin)
%ATFITTUNE fits linear tunes by scaling 2 quadrupole families
% NEWRING = ATFITTUNE(RING,NEWTUNES,QUADFAMILY1,QUADFAMILY2)
% NEWRING = ATFITTUNE(RING,DPP,NEWTUNES,QUADFAMILY1,QUADFAMILY2)
%
%RING:          Cell array
%DPP:           Optional momentum deviation (default 0)
%NEWTUNES:      Desired tune values
%QUADFAMILY1:   1st quadrupole family
%QUADFAMILY2:   2nd quadrupole family
%
%QUADFAMILY may be:
%   string: Family name
%   logical array: mask of selected elements in RING
%   Numeric array: list of selected elements in RING
%   Cell array: All elements selected by each cell
%
% NEWRING = ATFITTUNE(RING,...,'UseIntegerPart') With this flag, the 
% function fits the tunes to the total values of NEWTUNES, including 
% the integer part.
% With this option the function is substantially slower!
%

[UseIntegerPart,varargin]=getflag(varargin,'UseIntegerPart');
if isscalar(varargin{1}) && isnumeric(varargin{1})
    dpp=varargin{1};
    [newtunes,famname1,famname2]=deal(varargin{2:end});
else
    dpp=0;
    [newtunes,famname1,famname2]=deal(varargin{:});
end

idx1=varelem(ring,famname1);
idx2=varelem(ring,famname2);

kl1=atgetfieldvalues(ring(idx1),'PolynomB',{2});
kl2=atgetfieldvalues(ring(idx2),'PolynomB',{2});
delta = 1e-6;
if ~UseIntegerPart
    if any(newtunes>=1)
        warning('AT:FitTune','The integer part of the tunes is ignored unless you use the ''UseIntegerPart'' flag');
        newtunes=newtunes-floor(newtunes);
    end
    % Compute initial tunes before fitting
    [lindata, tunes] = atlinopt(ring,dpp); %#ok<ASGLU>
    % Take Derivative
    [lindata, tunes1] = atlinopt(setqp(ring,idx1,kl1,delta),dpp); %#ok<ASGLU>
    [lindata, tunes2] = atlinopt(setqp(ring,idx2,kl2,delta),dpp); %#ok<ASGLU>
else
    % Compute initial tunes before fitting
    lastpos=length(ring)+1;
    lindata = atlinopt(ring,dpp,1:lastpos);
    tunes=lindata(lastpos).mu/2/pi;
    % Take Derivative
    lindata1 = atlinopt(setqp(ring,idx1,kl1,delta),dpp,1:lastpos);
    lindata2 = atlinopt(setqp(ring,idx2,kl2,delta),dpp,1:lastpos);
    tunes1=lindata1(lastpos).mu/2/pi;
    tunes2=lindata2(lastpos).mu/2/pi;
end
%Construct the Jacobian
J = ([tunes1(:) tunes2(:)] - [tunes(:) tunes(:)])/delta;
dK = J\(newtunes(:)-tunes(:));

newring = setqp(ring,idx1,kl1,dK(1));
newring = setqp(newring,idx2,kl2,dK(2));

    function ring2=setqp(ring,idx,k0,delta)
        k=k0*(1+delta);
        ring2=ring;
        ring2(idx)=atsetfieldvalues(atsetfieldvalues(ring2(idx),'K',k),'PolynomB',{2},k);
    end

    function res=varelem(ring,arg)
        if islogical(arg)
            res=arg;
        elseif isnumeric(arg)
            res=false(size(ring));
            res(arg)=true;
        elseif ischar(arg)
            res=atgetcells(ring,'FamName',arg);
        elseif iscell(arg)
            res=false(size(ring));
            for i=1:length(arg)
                res=res|varelem(ring,arg{i});
            end
        else
            error('AT:GetElemList:WrongArg','Cannot parse argument');
        end
    end
end
