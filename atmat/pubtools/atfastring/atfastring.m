function [newring,newringrad] = atfastring(ring)
%ATFASTRING Generate simplified AT structures
%
%   [FASTRING,FASTRINGRAD]=ATFASTRING(RING)
%
%RING:          original AT structure, with no RF and no radiation.
%               Cavities are assumed to be at the beginning of the lattice
%
%FASTRING:      Strucuture containing a cavity, a linear 6x6 matrix and a
%               non-linear element simulating linear chromaticities and tune
%               shift with amplitudes
%
%FASTRINGRAD:   Strucuture containing a cavity, a diffusion element, a
%               damping matrix, a linear 6x6 transfer matrix and a
%               non-linear element simulating linear chromaticities and tune
%               shift with amplitudes

[~,nbper,voltage,harmnumber,U0]=atenergy(ring);
circumference=nbper*findspos(ring,length(ring)+1);
timelag=circumference*asin(U0/voltage)/2/pi/harmnumber;
[lindata,tunes,xsi]=atlinopt(ring,0);
gamma=(1+lindata.alpha.*lindata.alpha)./lindata.beta;

cavities=atgetcells(ring,'Frequency');
ignore=atgetcells(ring,'Class','RingParam');
frst=find(~(cavities | ignore),1);

figure(1);
orbit4=zeros(6,2);
orbit4(1:5,:)=findsyncorbit(ring,0,[1 frst]);
r1=detuning(ring,gamma,tunes,0.0005,orbit4(:,1));
%r1=detuning(ring,gamma,tunes,0.0001,orbit4(:,1));

figure(2);
%ringrad=atsetcavity(atradon(ring),6e6,true,992);
ringrad=atradon(ring);
orbit6=findorbit6(ringrad,[1 frst]);
r2=detuning(ringrad,gamma,tunes,0.0005,orbit6(:,1));
%r2=detuning(ringrad,gamma,tunes,0.0001,orbit6(:,1));

m66norad=findm66(ring(frst:end),[],orbit4(:,2));
m66norad=symplectify(m66norad);
m66rad=findm66(ringrad(frst:end),[],orbit6(:,2));

lin_elem=atM66('lin_elem',m66norad,'T1',-orbit4(:,2),'T2',orbit4(:,1),'Length',circumference);

damp_elem=atM66('damp_elem',m66norad\m66rad,'T1',-orbit6(:,2));
lin_elemrad=atM66('lin_elem',m66norad,'T2',orbit6(:,1),'Length',circumference);

quant_elem=atQuantDiff('quant_elem',ring);
nonlin_elem=atbaselem('nonlin_elem','DeltaQPass',...
    'Betax',lindata.beta(1),'Betay',lindata.beta(2),...
    'Alphax',lindata.alpha(1),'Alphay',lindata.alpha(2),...
    'Qpx',xsi(1),'Qpy',xsi(2),...
    'A1',r1(1),'A2',r1(2),'A3',r1(4));
nonlin_elemrad=atbaselem('nonlin_elem','DeltaQPass',...
    'Betax',lindata.beta(1),'Betay',lindata.beta(2),...
    'Alphax',lindata.alpha(1),'Alphay',lindata.alpha(2),...
    'Qpx',xsi(1),'Qpy',xsi(2),...
    'A1',r2(1),'A2',r2(2),'A3',r2(4));
newring=[ring(cavities);{lin_elem;nonlin_elem}];
newringrad=[atsetfieldvalues(ringrad(cavities),'TimeLag',timelag);...
    {quant_elem;damp_elem;lin_elemrad;nonlin_elemrad}];

    function r=detuning(ring,gamma,tunes,ampl,orbit)
        l=linspace(0,ampl.*ampl,10);
        subplot(2,1,1);
        [nuxx,nuzx]=atnuampl(ring,sqrt(l/gamma(1)),1,orbit);
        subplot(2,1,2);
        [nuxz,nuzz]=atnuampl(ring,sqrt(l/gamma(2)),3,orbit);
        if tunes(1) > 0.5
            nuxx=1-nuxx;
            nuxz=1-nuxz;
        end
        if tunes(2) > 0.5
            nuzx=1-nuzx;
            nuzz=1-nuzz;
        end
        %r=2*[l' ones(length(l),1)]\[nuxx' nuzx' nuxz' nuzz']
        r=2*([nuxx-nuxx(1);nuzx-nuzx(1);nuxz-nuxz(1);nuzz-nuzz(1)]*l')./(l(ones(4,1),:)*l');
    end
end
