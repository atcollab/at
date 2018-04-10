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
% NEWRING = ATFITTUNE(RING,...,'UseIntegerPart') Whith this flag, the 
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
if ~UseIntegerPart
    if newtunes(1)>1 || newtunes(2)>1
        warning('you passed also the integer part of the tunes, but it is ignored when you don''t use the flag ''UseIntegerPart''');
    end
    newtunes=newtunes-floor(newtunes);
end

kl1=atgetfieldvalues(ring(idx1),'PolynomB',{2});
kl2=atgetfieldvalues(ring(idx2),'PolynomB',{2});
if true
    delta = 1e-6;
    if ~UseIntegerPart
        % Compute initial tunes before fitting
        [lindata, tunes] = atlinopt(ring,dpp); %#ok<ASGLU>
        % Take Derivative
        [lindata, tunes1] = atlinopt(setqp(ring,idx1,kl1,delta),dpp); %#ok<ASGLU>
        [lindata, tunes2] = atlinopt(setqp(ring,idx2,kl2,delta),dpp); %#ok<ASGLU>
    else
        % Compute initial tunes before fitting
        lastpos=length(ring)+1;
        lindata = atlinopt(ring,dpp,1:lastpos); %#ok<ASGLU>
        tunes=lindata(lastpos).mu/2/pi;
        % Take Derivative
        lindata1 = atlinopt(setqp(ring,idx1,kl1,delta),dpp,1:lastpos); %#ok<ASGLU>
        lindata2 = atlinopt(setqp(ring,idx2,kl2,delta),dpp,1:lastpos); %#ok<ASGLU>
        tunes1=lindata1(lastpos).mu/2/pi;
        tunes2=lindata2(lastpos).mu/2/pi;
    end
    %Construct the Jacobian
    J = ([tunes1(:) tunes2(:)] - [tunes(:) tunes(:)])/delta;
    dK = J\(newtunes(:)-tunes(:));
else
    dK=0.01*fminsearch(@funtune,[0;0],...
        optimset(optimset('fminsearch'),'Display','iter','TolX',1.e-5));
end
newring = setqp(ring,idx1,kl1,dK(1));
newring = setqp(newring,idx2,kl2,dK(2));

%     function c=funtune(dK)
%         ring2=ring;
%         km1=kl1*(1+0.01*dK(1));
%         ring2(idx1)=atsetfieldvalues(atsetfieldvalues(ring2(idx1),'K',km1),'PolynomB',{2},km1);
%         km2=kl2*(1+0.01*dK(2));
%         ring2(idx2)=atsetfieldvalues(atsetfieldvalues(ring2(idx2),'K',km2),'PolynomB',{2},km2);
%         [lindata,tunes]=atlinopt(ring2,dpp); %#ok<SETNU>
%         dt=abs(newtunes(:)-tunes(:));
%         c=sum(dt.*dt);
%     end

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
