function r=PadPolynomAB(r)
% pads PolynomB and A to have the same number of elements


    r=cellfun(@(a)padpol(a),r,'un',0);

return

function a=padpol(a)

if isfield(a,'PolynomB')
try
    lpa=length(a.PolynomA);
    lpb=length(a.PolynomB);
catch
    a.PolynomA=[0];
    a.PolynomB=[0];
    a.MaxOrder=0;
    a.NumIntSteps=1;
    lpa=length(a.PolynomA);
    lpb=length(a.PolynomB);
end

    if lpa<lpb
    a.PolynomA=[a.PolynomA,zeros(1,lpb-lpa)];
    elseif lpa>lpb
    a.PolynomB=[a.PolynomB,zeros(1,lpa-lpb)];
    end
    
end

return