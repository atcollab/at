function varargout=plotB0curlyh(varargin)
%PLOTB0CURLYH  Plot B and H
%
% USAGE:
% >> atbaseplot(ring,@PLOTB0CURLYH);
% >> atplot(ring,@PLOTB0CURLYH);     (obsolete)
%
%  See also atbaseplot

if nargout == 1 % From atplot
    [lindata,ring,~]=deal(varargin{1:3});
    spl=PhysConstant.speed_of_light_in_vacuum.value;
    E0 = atenergy(ring);
    Brho=E0/spl;

    H=CurlyHlindata(lindata);
    B0=zeros(size(H));
    dips=findcells(ring,'BendingAngle');
    B0(dips)=Brho*getcellstruct(ring,'BendingAngle',dips)./...
        getcellstruct(ring,'Length',dips);


    plotdata(1).values=B0;
    plotdata(1).labels={['B [T] ' num2str(E0*1e-9) 'GeV']};
    plotdata(1).axislabel='B [T]';
    dispersion=cat(2,lindata.Dispersion)'; % right axis

    plotdata(2).values=[dispersion(:,1)*100 H*10000];
    plotdata(2).labels={'\eta_x [cm]','H [10-2]'};
    plotdata(2).axislabel='dispersion , H [cm]';
    varargout={plotdata};
else % From atbaseplot
    [ring,dpp]=deal(varargin{1:2});
    [linargs,varargs]=linoptions(varargin(3:end));
    [~,lindata]=atlinopt6(ring,1:length(ring)+1,linargs{:});
    s=cat(1,lindata.SPos);
    varargout={s,plotB0curlyh(lindata,ring,dpp,varargs{:})};
end
end
