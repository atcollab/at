function newring=atfitchrom(ring,varargin)
%ATFITCHROM Fit chromaticites by scaling 2 sextupole families
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
%NEWRING = ATFITTUNE(RING,...,'DPStep',dpstep)
%   dpstep is the momentum step applied to compute the chromaticity.
%   The default is the DPStep global option, which defaults to 3.0e-6
%
%NEWRING = ATFITTUNE(RING,...,'HStep',hstep)
%   hstep is the sextupole strength applied to build the jacobian [m^-3].
%   Default: 1.0e-5
%
% See also ATFITTUNE
allargs=getdparg(varargin);
newring=wrapper6d(ring,@xfit,allargs{:});

    function newring=xfit(ring,~,varargin)
        [deltaS,varargs]=getoption(varargin,'HStep',1.0e-5);
        [newchrom,famname1,famname2,varargs]=getargs(varargs,[],[],[]);

        idx1=varelem(ring,famname1);
        idx2=varelem(ring,famname2);
        kl1=atgetfieldvalues(ring(idx1),'PolynomB',{3});
        kl2=atgetfieldvalues(ring(idx2),'PolynomB',{3});
        %if true
        % Compute initial tunes before fitting
        chrom=getchrom(ring,varargs{:});

        % Take Derivative
        chrom1 = getchrom(setsx(ring,idx1,kl1,deltaS),varargs{:});
        chrom2 = getchrom(setsx(ring,idx2,kl2,deltaS),varargs{:});

        %Construct the Jacobian
        J = [chrom1-chrom chrom2-chrom]/deltaS;
        dK = J\(newchrom(:)-chrom);
        % else
        %     dK=fminsearch(@funchrom,[0;0],...
        %         optimset(optimset('fminsearch'),'Display','iter',...
        % 		'TolX',1.e-5,'TolFun',1.e-8));
        % end
        newring=setsx(ring,idx1,kl1,dK(1));
        newring=setsx(newring,idx2,kl2,dK(2));

        %     function c=funchrom(dK)
        %         ring2=ring;
        %         ring2(idx1)=atsetfieldvalues(ring2(idx1),'PolynomB',{3},kl1*(1+dK(1)));
        %         ring2(idx2)=atsetfieldvalues(ring2(idx2),'PolynomB',{3},kl2*(1+dK(2)));
        %         chrom = getchrom(ring2,dpp,deltaP);
        %         dt=abs(newchrom(:)-chrom(:));
        %         c=sum(dt.*dt);
        %     end

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
