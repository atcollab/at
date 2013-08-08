function newring=atfitchrom(ring,newchrom,famname1,famname2,varargin)
%ATFITTUNE fits chromaticites by scaling 2 sextupol families
% NEWRING = ATFITCHROM(RING,NEWCHROM,SEXTFAMILY1,SEXTFAMILY2)
%
%RING:          Cell array
%NEWCHROM:      Desired non-normalized chromaticities
%SEXTFAMILY1:   1st sextupole family
%SEXTFAMILY2:   2nd sextupole family
%
%SEXTFAMILY may be:
%   string: Family name
%   logical array: mask of selected elements in RING
%   Numeric array: list of selected elements in RING
%   Cell array: All elements selected by each cell

deltaP = 1e-8;
idx1=varelem(ring,famname1);
idx2=varelem(ring,famname2);
kl1=atgetfieldvalues(ring(idx1),'PolynomB',{3});
kl2=atgetfieldvalues(ring(idx2),'PolynomB',{3});
if true
    deltaS = 1e-5; % step size in Sextupole strngth
    
    % Compute initial tunes before fitting
    chrom=getchrom(ring,deltaP);
    
    % Take Derivative
    chrom1 = getchrom(setsx(ring,idx1,kl1,deltaS),deltaP);
    chrom2 = getchrom(setsx(ring,idx2,kl2,deltaS),deltaP);
    
    %Construct the Jacobian
    J = ([chrom1(:) chrom2(:)] - [chrom(:) chrom(:)])/deltaS;
    dK = J\(newchrom(:)-chrom(:));
else
    dK=fminsearch(@funchrom,[0;0],...
        optimset(optimset('fminsearch'),'Display','iter',...
		'TolX',1.e-5,'TolFun',1.e-8));
end
newring=setsx(ring,idx1,kl1,dK(1));
newring=setsx(newring,idx2,kl2,dK(2));

    function c=funchrom(dK)
        ring2=ring;
        ring2(idx1)=atsetfieldvalues(ring2(idx1),'PolynomB',{3},kl1*(1+dK(1)));
        ring2(idx2)=atsetfieldvalues(ring2(idx2),'PolynomB',{3},kl2*(1+dK(2)));
        chrom = getchrom(ring2,deltaP);
        dt=abs(newchrom(:)-chrom(:));
        c=sum(dt.*dt);
    end

    function ring2=setsx(ring,idx,k0,delta)
        ring2=atsetfieldvalues(ring,idx,'PolynomB',{3},k0*(1+delta));
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

    function chrom=getchrom(ring,dpp)
        [lindata, tunesa] = atlinopt(ring,0); %#ok<ASGLU>
        [lindata, tunesb] = atlinopt(ring,dpp); %#ok<ASGLU>
        chrom = (tunesb-tunesa)/deltaP;
    end
end
