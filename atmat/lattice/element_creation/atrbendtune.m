function elem=attunerbend(elem)
%ATTUNERBEND    Set X0ref and RefDZ for rectangular bending magnets
%
%NEWELEM=ATTUNERBEND(ELEM)
%   Set the X0ref and RefDZ attributes for rectangular bending magnets
%
%This function must be called after creating a rectangular bending magnet
%or after setting its polynomA/B attributes. It will set the correct X0ref
%and RefDZ attributes to get a zero closed orbit for the reference particle.
%
%Example:
%
%>> % Identify the rectangular bends
%>> rbends=atgetcells(ring,...);
%>> % Set their correct attributes
%>> ring(rbends)=cellfun(@attunerbend,ring(rbends),'UniformOutput',false);
%
%Does nothing if the passmethod is not a rectangular bend passmethod

% Remove radiation
passmethod=strrep(elem.PassMethod,'RadPass','Pass');

% Skip for non-rectangular magnets
if any(strcmp(passmethod,{'BndStrMPoleSymplectic4Pass', ...
        'ExactRectangularBendPass','ExactRectBendPass'}))
    theta=elem.BendingAngle;

    % Analytical estimate
    x0ref=elem.Length*((cos(0.5*theta)-1.0)/theta + sin(0.5*theta)/12);

    % Search if there are multipoles
    if checkmul(elem)
        x0ref=fzero(@cross,x0ref);
    end

    elem.X0ref=x0ref;
    vout=feval(passmethod,elem,zeros(6,1));
    elem.RefDZ=vout(6);
end

    function v=cross(x0r)
        % Return the horizontal exit angle of the reference particle
        elem.X0ref=x0r;
        rout=feval(passmethod,elem,zeros(6,1));
        v=rout(2);
    end

    function mult=checkmul(elem)
        % Check if there are multipoles
        mult=false;
        for order=1:elem.MaxOrder+1
            if elem.PolynomB(order) ~= 0.0
                mult=true;
                return;
            end
        end
    end

end