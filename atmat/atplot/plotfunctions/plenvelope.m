function varargout=plenvelope(varargin)
%PLENVELOPE    Plot beam envelope
%
%Helper function for atplot: plot
%- H and V beam envelopes on left axis

if nargout == 1
    refpts=1:length(varargin{2})+1;
    lind=atx(varargin{2:3},refpts);
    sigma=cat(3,lind.beam66);
    plotdata(1).values=1.0e6*sqrt([squeeze(sigma(1,1,:)) squeeze(sigma(3,3,:))]);
    plotdata(1).labels={'horizontal','vertical'};
    plotdata(1).axislabel='beam envelope [\mum]';
    varargout={plotdata};
else
    s=findspos(varargin{1},1:length(varargin{1})+1);
    varargout={s,plenvelope([],varargin{1:2})};
end
end
