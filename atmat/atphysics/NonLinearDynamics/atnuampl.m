function varargout=atnuampl(ring,ampl,xz,varargin)
%ATNUAMPL	computes tune shift with amplitude
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE)
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE,1)
%
%	Computes tunes for the specified horizontal amplitudes
%
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE,3)
%
%	Computes tunes for the specified vertical amplitudes
%
%ATNUAMPL(...)
%   Plots the computed tunes in the current axes
%
%ATNUAMPL(...,Name,Value)
%   Uses additional options specified by one or more Name,Value pairs.
%   Possible values are:
%       nturns: specify the number of turns for tracking (default 256)
%       method: specify the method for tune determination
%               1: Highest peak in fft
%               2: Interpolation on fft results
%               3: Windowing + interpolation (default)
%               4: NAFF
%   Other options are transmitted to the plot function


lab={'x^2','p_x^2','z^2','p_z^2'};
if nargin < 3, xz=1; end
if ~isempty(varargin) && isnumeric(varargin{1})
    orbit=varargin{1};
    varargin(1)=[];
else
    warning off MATLAB:singularMatrix
    orbit=findorbit6(ring);
    warning on MATLAB:singularMatrix
    if ~all(isfinite(orbit))
        orbit=zeros(6,1);
        orbit(1:5)=findsyncorbit(ring,0);
    end
end
[nturns,varargin]=getoption(varargin,'nturns',256);
[method,varargin]=getoption(varargin,'method',3);

[~,nbper]=atenergy(ring);
[lindata,fractune0]=atlinopt(ring,0,1:length(ring)+1);
tune0=nbper*lindata(end).mu/2/pi;
offs=[nbper -nbper];
siza=size(ampl);
nampl=prod(siza);
p0=repmat(0.00003*[1;0;1;0;0;0], 1,nampl); % 30 microns minimum amplitude
p0(xz,:)=max(p0(xz,:),ampl(:)');
p0=p0+orbit(:,ones(1,nampl));
p1=ringpass(ring,p0,nturns)-orbit(:,ones(1,nampl*nturns));
tunetrack=[findtune(reshape(p1(1,:),nampl,nturns)',method,reshape(p1(2,:),nampl,nturns)');...
    findtune(reshape(p1(3,:),nampl,nturns)',method,reshape(p1(4,:),nampl,nturns)')]';
[~,k]=min([fractune0-tunetrack(1,:); 1-fractune0-tunetrack(1,:)]);
np=offs(k);
offset=round(tune0-np.*tunetrack(1,:));
tunetrack=np(ones(nampl,1),:).*tunetrack + offset(ones(nampl,1),:);
if nargout > 0
    varargout={reshape(tunetrack(:,1),siza),reshape(tunetrack(:,2),siza)};
else
    inttunes=floor(tune0);
    plot((ampl.*ampl)',tunetrack-inttunes(ones(nampl,1),:),'o-',varargin{:});
    legend('\nu_x','\nu_z');
    xlabel(lab{xz});
    ylabel('\nu');
    grid on
end
