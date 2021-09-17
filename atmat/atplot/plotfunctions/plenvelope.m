function varargout=plenvelope(varargin)
%PLENVELOPE    Plot beam envelope
%
%Helper function for atplot: plot
%- H and V beam envelopes on left axis
%
% USAGE:
% >> atbaseplot(ring,@plenvelope);
% >> atplot(ring,@plenvelope);     (obsolete)
%
%  See also atbaseplot

if nargout == 1 % From atplot
    [~,ring,dpp]=deal(varargin{1:3});
    lind=atx(ring,dpp,1:length(ring)+1);
    sigma=cat(3,lind.beam66);
    plotdata(1).values=1.0e6*sqrt([squeeze(sigma(1,1,:)) squeeze(sigma(3,3,:))]);
    plotdata(1).labels={'horizontal','vertical'};
    plotdata(1).axislabel='beam envelope [\mum]';
    varargout={plotdata};
else % From atbaseplot
    s=findspos(varargin{1},1:length(varargin{1})+1);
    varargout={s,plenvelope([],varargin{:})};
end
end
