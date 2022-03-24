function varargout=plCorrectorStrength(varargin)
%PLCORRECTORSTRENGTH Plot PolynomB
%Helper function for atplot: plot
%- PolynomB(1), PolynomA(1), PolynomA(2), PolynomA(2) on the left axis
%
%  EXAMPLEs
% >> atbaseplot(ring,'synopt',false,@plCorrectorStrength);
% >> atplot(ring,@plCorrectorStrength,'synopt',false);     (obsolete)
%
%  See also atplot atbaseplot
%

if nargout == 1 % From atplot
    [~,ring]=deal(varargin{1:2});
    sz=length(ring)+1;

    PolynomBVal1=zeros(sz,1);
    PolynomBVal2=zeros(sz,1);
    PolynomBVal3=zeros(sz,1);
    PolynomBVal4=zeros(sz,1);
    ind=findcells(ring,'PolynomB');

    PolynomBVal1(ind)=getcellstruct(ring,'PolynomB',ind,1,1);
    PolynomBVal2(ind)=getcellstruct(ring,'PolynomA',ind,1,1);
    PolynomBVal3(ind)=getcellstruct(ring,'PolynomA',ind,1,2);
    PolynomBVal4(ind)=getcellstruct(ring,'PolynomB',ind,1,2);


    plotdata(1).values=[PolynomBVal1 PolynomBVal2 ...
        PolynomBVal3/1000 PolynomBVal4/1000 ];
    plotdata(1).labels={'H cor','V cor',...
        'Skew Quad/1000','Quad/1000'};
    plotdata(1).axislabel='PolynomB';
    varargout={plotdata};
else % From atbaseplot
    [ring,dpp]=deal(varargin{1:2});
    s=findspos(ring,1:length(ring)+1)';
    varargout={s,plCorrectorStrength([],ring,dpp,varargin{3:end})};
end

end
