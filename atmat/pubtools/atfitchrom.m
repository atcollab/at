function newring=atfitchrom(ring,newchrom,famname1,famname2,varargin)
%ATFITTUNE fits chromaticites by scaling 2 sextupol families
% NEWRING = ATFITCHROM(RING,NEWCHROM,SEXTFAMILY1,SEXTFAMILY2)

deltaP = 1e-8;
idx1=findfam(famname1);
idx2=findfam(famname2);
kl1=getcellstruct(ring,'PolynomB',idx1,3);
kl2=getcellstruct(ring,'PolynomB',idx2,3);
if true
    deltaS = 1e-5; % step size in Sextupole strngth
    
    % Compute initial tunes before fitting
    chrom=getchrom(ring,deltaP);
    
    % Take Derivative
    newring=setcellstruct(ring,'PolynomB',idx1,kl1*(1+deltaS),3);
    chrom1 = getchrom(newring,deltaP);
    newring=setcellstruct(ring,'PolynomB',idx2,kl2*(1+deltaS),3);
    chrom2 = getchrom(newring,deltaP);
    
    %Construct the Jacobian
    J = ([chrom1(:) chrom2(:)] - [chrom(:) chrom(:)])/deltaS;
    dK = J\(newchrom(:)-chrom(:));
else
    dK=fminsearch(@funchrom,[0;0],...
        optimset(optimset('fminsearch'),'Display','iter',...
		'TolX',1.e-5,'TolFun',1.e-8));
end
newring=setcellstruct(ring ,'PolynomB',idx1,kl1*(1+dK(1)),3);
newring=setcellstruct(newring,'PolynomB',idx2,kl2*(1+dK(2)),3);

    function c=funchrom(dK)
        ring2=setcellstruct(ring ,'PolynomB',idx1,kl1*(1+dK(1)),3);
        ring2=setcellstruct(ring2,'PolynomB',idx2,kl2*(1+dK(2)),3);
        chrom = getchrom(ring2,deltaP);
        dt=abs(newchrom(:)-chrom(:));
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

    function chrom=getchrom(ring,dpp)
        [lindata, tunesa] = linopt(ring,0);
        [lindata, tunesb] = linopt(ring,dpp);
        chrom = (tunesb-tunesa)/deltaP;
    end
end
