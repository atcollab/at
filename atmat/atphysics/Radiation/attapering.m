function ring = attapering(ring,varargin)
%ATTAPERING Scale magnet strengths
%
%NEWRING=ATTAPERING(RING)   Scales dipole strengths with local energy to
%   cancel the closed orbit due to synchrotron radiation.
%
%NEWRING=ATTAPERING(RING,'multipoles')  Scales also the multipoles to
%   cancel optics errors
%
%NEWRING=ATTAPERING(RING,'niter',niter) Performs niter iterations (useful
%   when multipoles are scaled)
%

[niter, varargs] = getoption(varargin, 'niter', 1);
[multipoles, varargs] = getflag(varargs, 'multipoles'); %#ok<ASGLU>
dipin = atgetcells(ring, 'BendingAngle'); % Dipole entrance
dipout = circshift(dipin, 1);             % Dipole exit
    multin = atgetcells(ring,'PolynomB') & ~dipin;
    multout = circshift(multin, 1);
for it=1:niter
    o6 = findorbit6(ring,1:length(ring)+1);
    dppd = 0.5*(o6(5, dipin) + o6(5, dipout));
    ring = atsetfieldvalues(ring,dipin,'FieldScaling',1.0+dppd');
    
    if multipoles
        dppm = 0.5*(o6(5, multin) + o6(5, multout));
        ring(multin) = cellfun(@scale,ring(multin),num2cell(dppm'), ...
            'UniformOutput', false);
    end
end

    function elem=scale(elem, factor)
        elem.PolynomB=factor*elem.PolynomB;
        elem.PolynomQ=factor*elem.PolynomA;
    end

end