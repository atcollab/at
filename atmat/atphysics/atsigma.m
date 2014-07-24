function sigma = atsigma(varargin)
%ATSIGMA Construct a abeam sigma matrix
%
%   SIGMA=ATSIGMA(BETA,ALPHA,EMIT)
%       builds a 2x2 sigma matrix for a transverse plane
%
%   SIGMA=ATSIGMA(ESPREAD,BLENGTH)
%       builds a 2x2 sigma matrix for the longitudinal plane
%
%   SIGMA=ATSIGMA(BETAX,ALPHAX,EMITX,BETAZ,ALPHAZ,EMITZ)
%       builds a 4x4 sigma matrix
%
%   SIGMA=ATSIGMA(BETAX,ALPHAX,EMITX,BETAZ,ALPHAZ,EMITZ,ESPREAD,BLENGTH)
%       builds a 6x6 sigma matrix
%
%   SIGMA=ATSIGMA(ATSTRUCT)
%       builds a 6x6 sigma matrix

if iscell(varargin{1})
    beamdata=atx(varargin{1},0,1);
    sigma=beamdata(1).beam66;
elseif nargin == 2
    [espread,blength]=deal(varargin{:});
    sigma=[espread.*espread 0;0 blength.*blength];
elseif nargin == 3
    [bx,ax,epsx]=deal(varargin{:});
    sigma=epsx*[bx -ax;-ax (1+ax*ax)/bx];
elseif nargin == 6
    [bx,ax,epsx,bz,az,epsz]=deal(varargin{:});
    sigma=[atsigma(bx,ax,epsx) zeros(2,2);...
        zeros(2,2) atsigma(bz,az,epsz)];
else
    [bx,ax,epsx,bz,az,epsz,espread,blength]=deal(varargin{:});
    sigma=[atsigma(bx,ax,epsx) zeros(2,4);...
        zeros(2,2) atsigma(bz,az,epsz) zeros(2,2);...
        zeros(2,4) atsigma(espread,blength)];
end
end
