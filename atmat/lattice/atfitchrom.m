function newring=atfitchrom(ring,varargin)
%ATFITCHROM Fit chromaticites by scaling 2 sextupole families
%
% NEWRING = ATFITCHROM(RING,NEWCHROM,SEXTFAMILY1,SEXTFAMILY2)
% NEWRING = ATFITCHROM(RING,DPP,NEWCHROM,SEXTFAMILY1,SEXTFAMILY2)
%
%RING:          Cell array
%DPP:           Optional momentum deviation (default 0)
%NEWCHROM:      Desired non-normalized chromaticities
%SEXTFAMILY1:   1st sextupole family
%SEXTFAMILY2:   2nd sextupole family
%
%SEXTFAMILY may be:
%   string: Family name
%   logical array: mask of selected elements in RING
%   Numeric array: list of selected elements in RING
%   Cell array: All elements selected by each cell
%
%NEWRING = ATFITCHROM(RING,...,'DPStep',dpstep)
%   dpstep is the momentum step applied to compute the chromaticity.
%   The default is the DPStep global option, which defaults to 3.0e-6
%
%NEWRING = ATFITCHROM(RING,...,'HStep',hstep)
%   hstep is the sextupole strength applied to build the jacobian [m^-3].
%   Default: 1.0e-5
%
%NEWRING = ATFITCHROM(RING,...,'dp',DP)
%   Specify off-momentum if radiation is OFF (default 0)
%
%NEWRING = ATFITCHROM(RING,...,'dct',DCT)
%   Specify path lengthening if radiation is OFF (default 0)
%
% See also ATFITTUNE

newargs=getdparg(varargin);
newring=wrapper6d(ring,@fitchrom,newargs{:});

    function newring=fitchrom(ring,~,varargin)
        [deltaS,varargs]=getoption(varargin,'HStep',1.0e-5);
        [newchrom,famname1,famname2,varargs]=getargs(varargs,[],[],[]);

        idx1=varelem(ring,famname1);
        idx2=varelem(ring,famname2);
        kl1=atgetfieldvalues(ring(idx1),'PolynomB',{3});
        kl2=atgetfieldvalues(ring(idx2),'PolynomB',{3});
        %if true
        % Compute initial chromaticities before fitting
        chrom=getchrom(ring,varargs{:});

        % Take Derivative
        chrom1 = getchrom(setsx(ring,idx1,kl1,deltaS),varargs{:});
        chrom2 = getchrom(setsx(ring,idx2,kl2,deltaS),varargs{:});

        % Construct the Jacobian
        J = [chrom1-chrom chrom2-chrom]/deltaS;
        dK = J\(newchrom(:)-chrom);

        % Apply new strengths
        newring=setsx(ring,idx1,kl1,dK(1));
        newring=setsx(newring,idx2,kl2,dK(2));

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

        function chrom=getchrom(ring,varargin)
            [~,chr]=tunechrom(ring,varargin{:});
            chrom=chr(1:2)';
        end
    end
end
