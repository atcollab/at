function varargout=plotAperture(varargin) 
%PLOTAPERTURE Plots x and y aperture.
%
%Helper function for atplot: plot the physical aperture
%
%  USAGE:
% >> atbaseplot(ring,@plotAperture);
% >> atplot(ring,@plotAperture);        (obsolete)
%
%See also atplot atbaseplot

if nargout == 1 % From atplot
    ring=varargin{2};
    apind=findcells(ring,'Limits');
    
    xm=getcellstruct(ring,'Limits',apind,1);
    ym=getcellstruct(ring,'Limits',apind,3);
    xp=getcellstruct(ring,'Limits',apind,2);
    yp=getcellstruct(ring,'Limits',apind,4);
    
    Xp=[nan(size(ring)); nan];
    Xm=Xp;
    Yp=Xp;
    Ym=Xp;
    
    Xp(apind)=xp;
    Xm(apind)=xm;
    Yp(apind)=yp;
    Ym(apind)=ym;
    
    plotdata(1).values=[Xp Xm Yp Ym]*1e2;%
    plotdata(1).labels={'x aperture','x aperture','y aperture','y aperture'};
    plotdata(1).axislabel='x or y aperture [cm]';
    %
    % plotdata(2).values=[Yp Ym]*1e2;%
    % plotdata(2).labels={'y aperture','y aperture'};
    % plotdata(2).axislabel='y aperture [cm]';
    varargout={plotdata};
else                % From atbaseplot
    s=findspos(varargin{1},1:length(varargin{1})+1);
    varargout={s,plotAperture([],varargin{:})};
end

end 