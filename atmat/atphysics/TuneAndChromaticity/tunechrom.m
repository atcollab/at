function varargout = tunechrom(ring,varargin)
%TUNECHROM computes linear tunes and chromaticities
%
%TUNE=TUNECHROM(RING)	Quick calculation of the fractional part of the tune
%	from numerically computed transfer matrix.
%
% [TUNE, CHROM] = TUNECHROM(RINGD,DP,'get_chrom') - optionally computes the
%    chromaticities by numerical differentiation from the difference between
%   tune values at momentums DP+0.5*DPStep and DP-0.5*DPStep
%
%[...]=TUNECHROM(...,'orbit',ORBITIN	Do not search for closed orbit.
%   Instead ORBITIN,a 6x1 vector of initial conditions is used:
%   This syntax is useful to avoid recomputing the closed orbit if is
%   already known;
%
%[...]=TUNECHROM(RING,DP)       (obsolete)
%[...]=TUNECHROM(RING,...,'dp',DP)	Specify the momentum deviation when
%   radiation is OFF (default: 0)
%
%[...]=TUNECHROM(RING,...,'dct',DCT) Specify the path lengthening when
%   radiation is OFF (default: 0)
%
%[...]=TUNECHROM(RING,...,'df',DF) Specify the RF frequency deviation when
%   radiation is OFF (default: 0)
%
% Note: TUNECHROM computes tunes and chromaticities from the one-turn
%   transfer matrix. The transfer matrix is computed from tracking using
%   numerical differentiation. The error of numerical differentiation
%   is sensitive to the step size. (Reference: Numerical Recipes)
%   The calculation of tunes involves one numerical differentiation.
%   The calculation of chromaticity involves TWO!!! numerical differentiations.
%   The error in calculated chromaticity from may be substantial (~ 1e-5).
%   Use the XYStep and DPStep keyword arguments to control the step size
%   in chromaticity calculations
%
% See also ATLINOPT6
[cpl1,allargs]=getflag(varargin,'coupling');	% ignored, kept for compatibility
[cpl2,allargs]=getoption(allargs,'coupled',true);	% ignored, kept for compatibility
allargs=getdparg(allargs);
if cpl1 || ~cpl2
    warning('AT:ObsoleteParameter','The "coupled" flag is ignored: coupling is always assumed');
end
[varargout{1:nargout}]=frequency_control(@xtunechrom,ring,allargs{:});

    function [tune, chrom] = xtunechrom(ring,varargin)
        [oldchrom,varargs]=getflag(varargin,'chrom');
        [chrom,varargs]=getflag(varargs,'get_chrom');
        [dpargs,varargs]=getoption(varargs,{'orbit','dp','dct','df'});
        [cavargs,varargs]=getoption(varargs,{'cavpts'});
        [DPStep,~]=getoption(varargs,'DPStep');
        is_6d=getoption(varargs,'is_6d',[]); % Always set by frequency_control, keep in varargs
        if is_6d
            tunefunc=@tune6;
        else
            tunefunc=@tune4;
        end

        [~,orbitin]=findorbit(ring,dpargs{:},varargs{:});
        
        tune=tunefunc(ring,'orbit',orbitin,varargs{:});

        if chrom || oldchrom || nargout == 2
            if is_6d
                frf=get_rf_frequency(ring);
                DFStep=-DPStep*mcf(atradoff(ring))*frf;
                rgup=atsetcavity(ring,'Frequency',frf+0.5*DFStep,cavargs{:});
                rgdn=atsetcavity(ring,'Frequency',frf-0.5*DFStep,cavargs{:});
                [~,o1P]=findorbit6(rgup,'guess',orbitin,varargs{:});
                [~,o1M]=findorbit6(rgdn,'guess',orbitin,varargs{:});
                deltap=o1P(5)-o1M(5);
            else
                dp=orbitin(5);
                [~,o1P]=findorbit4(ring,dp+0.5*DPStep,'guess',orbitin,varargs{:},'strict',-1);
                [~,o1M]=findorbit4(ring,dp-0.5*DPStep,'guess',orbitin,varargs{:},'strict',-1);
                deltap=DPStep;
            end
            tuneP=tunefunc(ring,'orbit',o1P,varargs{:});
            tuneM=tunefunc(ring,'orbit',o1M,varargs{:});
            chrom = (tuneP - tuneM)/deltap;
        end

        function f=get_rf_frequency(ring)
            % Get the initial RF frequency
            cavities=ring(atgetcells(ring, 'Frequency'));
            freqs=atgetfieldvalues(cavities,'Frequency');
            f=freqs(1);
        end

    function tune=tune4(ring,varargin)
        mm=findm44(ring,NaN,varargin{:});
        [~,vals]=amat(mm);
        tune=mod(angle(vals)/2/pi,1);
    end

    function tune=tune6(ring,varargin)
        mm=findm66(ring,varargin{:});
        [~,vals]=amat(mm);
        tune=mod(angle(vals)/2/pi,1);
    end
    end
end
