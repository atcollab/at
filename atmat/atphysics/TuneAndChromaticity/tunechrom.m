function [tune, varargout] = tunechrom(ring,varargin)
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
%[...]=TUNECHROM(...,'dp',DP)	Specify the momentum deviation when
%   radiation is OFF (default: 0)
%
%[...]=TUNECHROM(...,'dct',DCT) Specify the path lengthening when
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

[oldchrom,varargs]=getflag(varargin,'chrom');   
[chrom,varargs]=getflag(varargs,'get_chrom');
[cpl1,varargs]=getflag(varargs,'coupling');	% ignored, kept for compatibility
[cpl2,varargs]=getoption(varargs,'coupled',true);	% ignored, kept for compatibility
[orbitin,varargs]=getoption(varargs,'orbit',[]);
[dp,varargs]=getoption(varargs,'dp',NaN);
[dp,varargs]=getargs(varargs,dp,'check',@(arg) isscalar(arg) || isnumeric(arg));
[DPStep,~]=getoption(varargs,'DPStep');

if cpl1 || ~cpl2
    warning('AT:ObsoleteParameter','The "coupled" flag is ignored: coupling is always assumed');
end

[~,orbitin]=findorbit(ring,'dp',dp,'orbit',orbitin,varargs{:});
is6d=check_radiation(ring);
if is6d
    mm=findm66(ring,'orbit',orbitin,varargs{:});
else
    dp=orbitin(5);
    mm=findm44(ring,dp,'orbit',orbitin,varargs{:});
end
[~,vals]=amat(mm);
tune=mod(angle(vals)/2/pi,1);

if chrom || oldchrom
    if is6d
        frf=get_rf_frequency(ring);
        DFStep=-DPStep*mcf(atradoff(ring))*frf;
        rgup=atsetcavity(ring,'Frequency',frf+0.5*DFStep);
        rgdn=atsetcavity(ring,'Frequency',frf-0.5*DFStep);
        [~,o1P]=findorbit(rgup,[],'guess',orbitin,varargs{:});
        [~,o1M]=findorbit(rgdn,[],'guess',orbitin,varargs{:});
        tuneP=tunechrom(rgup,'orbit',o1P,varargs{:});
        tuneM=tunechrom(rgdn,'orbit',o1M,varargs{:});
        deltap=o1P(5)-o1M(5);
    else
        dp=orbitin(5);
        tuneP=tunechrom(ring,'dp',dp + 0.5*DPStep,varargs{:});
        tuneM=tunechrom(ring,'dp',dp - 0.5*DPStep,varargs{:});
        deltap=DPStep;
    end
    varargout{1} = (tuneP - tuneM)/deltap;
end
    
    function f=get_rf_frequency(ring)
        % Get the initial RF frequency
        cavities=ring(atgetcells(ring, 'Frequency'));
        freqs=atgetfieldvalues(cavities,'Frequency');
        f=freqs(1);
    end
end
