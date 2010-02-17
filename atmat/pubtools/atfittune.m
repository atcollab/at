function newring=atfittune(ring,newtunes,famname1,famname2,varargin)
%ATFITTUNE fits linear tunes by scaling 2 quadrupole families
% NEWRING = ATFITTUNE(RING,NEWTUNES,QUADFAMILY1,QUADFAMILY2)

idx1=findfam(famname1);
idx2=findfam(famname2);
kl1=getcellstruct(ring,'K',idx1);
kl2=getcellstruct(ring,'K',idx2);
if true
    delta = 1e-6;

    % Compute initial tunes before fitting
    [lindata, tunes] = linopt(ring,0);

    % Take Derivative
    [lindata, tunes1] = linopt(setqp(ring,idx1,kl1,delta),0);
    [lindata, tunes2] = linopt(setqp(ring,idx2,kl2,delta),0);

    %Construct the Jacobian
    J = ([tunes1(:) tunes2(:)] - [tunes(:) tunes(:)])/delta;
    dK = J\(newtunes(:)-tunes(:));
else
    dK=0.01*fminsearch(@funtune,[0;0],...
        optimset(optimset('fminsearch'),'Display','iter','TolX',1.e-5));
end
newring = setqp(ring,idx1,kl1,dK(1));
newring = setqp(newring,idx2,kl2,dK(2));

    function c=funtune(dK)
        ring2 = setqp(ring ,idx1,kl1,0.01*dK(1));
        ring2 = setqp(ring2,idx2,kl2,0.01*dK(2));
        [lindata,tunes]=linopt(ring2,0); %#ok<SETNU>
        dt=abs(newtunes(:)-tunes(:));
        c=sum(dt.*dt);
    end

    function idx=findfam(famname)
        if iscell(famname)
            idx=[];
            for i=1:numel(famname)
                idx=[idx findfam(famname{i})]; %#ok<AGROW>
            end
        else
            idx=findcells(ring,'FamName',famname);
        end
    end

    function ring2=setqp(ring,idx,k0,delta)
        k=k0*(1+delta);
        ring2=setcellstruct(ring,'K',idx,k);
        ring2=setcellstruct(ring2,'PolynomB',idx,k,2);
    end
end
