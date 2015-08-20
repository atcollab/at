function [newring,newringrad] = atfastring(ring0,split)
%ATFASTRING Generate simplified AT structures
%
% The given ring structure is modified so that cavities are moved to the
% beginning and the rest of the ring is replaced by a linear 6x6 transfer
% matrix followed by a non-linear element providing tune shifts with
% amplitude and momentum.
%
%   [FASTRING,FASTRINGRAD]=ATFASTRING(RING)
%
%RING:          original AT structure, with no RF and no radiation.
%               Cavities are assumed to be at the beginning of the lattice
%
%FASTRING:      Structure containing a cavity, a linear 6x6 matrix and a
%               non-linear element simulating linear chromaticities and
%               tune shift with amplitudes
%
%FASTRINGRAD:   Strucuture containing a cavity, a diffusion element,
%               a linear 6x6 transfer matrix and a non-linear element
%               simulating linear chromaticities and tune shift with
%               amplitudes
%
%   [FASTRING,FASTRINGRAD]=ATFASTRING(RING,REFPTS)
%
% The ring is split at the specified locations, and each section is
% transformed in the same way as previously described


global GLOBVAL

if nargin < 2
    split=[];
end
if islogical(split)
    iend=[find(split) length(ring0)+1];
else
    iend=[split length(ring0)+1];
end
ibeg=[1 iend(1:end-1)];

xm=0.001;
zm=0.0005;
GLOBVAL.E0=atenergy(ring0);
[lindata,tunes,xsi]=atlinopt(ring0,0); %#ok<ASGLU>
gamma=(1+lindata.alpha.*lindata.alpha)./lindata.beta;

ringv=arrayfun(@rearrange,ibeg,iend,'UniformOutput',false);
ring=cat(1,ringv{:});
markers=atgetcells(ring,'FamName','xbeg|xend');
ringrad=atradon(ring);

orbit4=zeros(6,sum(markers));
orbit4(1:5,:)=findsyncorbit(ring,0,markers);
orbit4=num2cell(orbit4,1);
r1=detuning(ring,gamma,xm,zm,orbit4(:,1));

orbit6=num2cell(findorbit6(ringrad,markers),1);

counter=0;
[rv,rvrad]=cellfun(@rebuild,ringv,orbit4(1:2:end),orbit6(1:2:end),...
    orbit4(2:2:end),orbit6(2:2:end),'UniformOutput',false);
nonlin_elem=atbaselem('NonLinear','DeltaQPass',...
    'Betax',lindata.beta(1),'Betay',lindata.beta(2),...
    'Alphax',lindata.alpha(1),'Alphay',lindata.alpha(2),...
    'Qpx',xsi(1),'Qpy',xsi(2),...
    'A1',r1(1),'A2',r1(2),'A3',r1(4),...
    'T1',-orbit4{end},'T2',orbit4{end});
nonlin_elemrad=atbaselem('NonLinear','DeltaQPass',...
    'Betax',lindata.beta(1),'Betay',lindata.beta(2),...
    'Alphax',lindata.alpha(1),'Alphay',lindata.alpha(2),...
    'Qpx',xsi(1),'Qpy',xsi(2),...
    'A1',r1(1),'A2',r1(2),'A3',r1(4),...
    'T1',-orbit6{end},'T2',orbit6{end});
newring=cat(1,rv{:},nonlin_elem);
newringrad=cat(1,rvrad{:},nonlin_elemrad);

    function rg=rearrange(i1,i2)
        slice=ring0(i1:i2-1);
        cav=atgetcells(slice,'Frequency') | atgetcells(slice,'Class','RingParam');
        rg=[slice(cav);atmarker('xbeg');slice(~cav);atmarker('xend')];
    end
    function [rg,rgrad]=rebuild(slice,o4b,o6b,o4e,o6e)
        counter=counter+1;
        cc=num2str(counter);
%         m1=atmarker(['xbeg' cc]);
%         m2=atmarker(['xend' cc]);
        i1=find(atgetcells(slice,'FamName','xbeg'),1);
        dipoles=atgetcells(slice,'BendingAngle');
        theta=atgetfieldvalues(slice(dipoles),'BendingAngle');
        lendp=atgetfieldvalues(slice(dipoles),'Length');
        s=diff(findspos(slice,[1 length(slice)+1]));
        I2=sum(abs(theta.*theta./lendp));
        
        m66norad=symplectify(findm66(slice(i1:end),[],o4b));
        lin_elem=atM66(['Linear_' cc],m66norad,'T1',-o4b,'T2',o4e,'Length',s,'I2',I2);
        rg=[slice(1:i1-1);lin_elem];
        
        [slicerad,radindex]=atradon(slice);
        diff_elem=atQuantDiff(['Diffusion_' cc],quantumDiff(slicerad,radindex,o6b));
        m66rad=findm66(slicerad(i1:end),[],o6b);
        lin_elemrad=atM66(['Linear_' cc],m66rad,'T1',-o6b,'T2',o6e,'Length',s,'I2',I2);
        rgrad=[slicerad(1:i1-1);diff_elem;lin_elemrad];
    end
    function r=detuning(ring,gamma,xm,zm,orbit)
        x2=linspace(0,xm.*xm,10);
        z2=linspace(0,zm.*zm,10);
        [nuxx,nuzx]=atnuampl(ring,sqrt(x2),1,orbit);
        [nuxz,nuzz]=atnuampl(ring,sqrt(z2),3,orbit);
        tune0=floor([nuxx(1);nuzz(1)]);
        subplot(2,1,1);
        plot(x2,[nuxx;nuzx]-tune0(:,ones(1,10)));
        subplot(2,1,2);
        plot(z2,[nuxz;nuzz]-tune0(:,ones(1,10)));
        rx=([nuxx-nuxx(1);nuzx-nuzx(1)]*x2')./(x2([1 1],:)*x2')/gamma(1);
        rz=([nuxz-nuxz(1);nuzz-nuzz(1)]*z2')./(z2([1 1],:)*z2')/gamma(2);
        r=2*[rx;rz];
    end
end
