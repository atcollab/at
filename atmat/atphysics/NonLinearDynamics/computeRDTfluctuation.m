function [RDT,buildup_fluctuation,natural_fluctuation] = computeRDTfluctuation(ring, varargin)
%COMPUTERDTFLUCTUATION Computes Hamiltonian resonance driving terms (RDTs)
%   This function is based on simplestoragering and returns the RDTs
%   and their longitudinal fluctuations.
%   
% [RDT,buildup_fluctuation,natural_fluctuation]=computeRDTfluctuation(ring, varargin)
%
%   ring is the AT lattice
%  The additional arguments:
%   nslices: number of slices of each sextupole, which affects the crossing
%       terms. default: 4.
%   nperiods: number of periods. RDTs and RDT build-up fluctuations will be
%       computed for n periods.  default: 1.
%       natural RDT fluctuation of different periods are the same.
%       So the results contain only one period.
%   
%   RDT: struct, RDTs (complex numbers) and 
%       amplitude-dependent tune shifts (real)
%       (ADTS are calculated using h22000, h11110 and h00220)
%   buildup_fluctuation: a struct of complex arrays, 
%       accumulated RDTs from s=0,
%       showing the build-up and cancellation of RDTs along the 
%       longitudinal position.
%   natural_fluctuation: a struct of double arrays,
%       absolute values of one-period RDTs observed at different
%       longitudinal starting position.
%       same as s_dependent_driving_terms in ELEGANT.
%   
%  REFERENCES
%    [1] Johan Bengtsson, SLS Note 09/97, (1997)
%    [2] S. C. Leemann, A. Streun, Phys. Rev. ST Accel. Beams 14, 030701 (2011)
%    [3] A. Franchi, L. Farvacque, F. Ewald, G. Le Bec, and K. B. Scheidt, Phys. Rev. ST Accel. Beams 17, 074001 (2014)
%    [4] B. Wei, Z. Bai, J. Tan, L. Wang, and G. Feng, Phys. Rev. Accel. Beams 26, 084001 (2023)
%

% Validates the input arguements.
[nslices,varargs]=getoption(varargin,'nslices',4);
[nperiods,~]=getoption(varargs,'nperiods',1);

if ~(isnumeric(nslices) && isreal(nslices) && isfinite(nslices) ...
        && (nslices > 0) && round(nslices) == nslices)
    throw(MException('RDTFluctuation:variableTypeError', ...
        'nslices must be a positive integer'))
end
if ~(isnumeric(nperiods) && isreal(nperiods) && isfinite(nperiods) ...
        && (nperiods > 0) && round(nperiods) == nperiods)
    throw(MException('RDTFluctuation:variableTypeError', ...
        'nperiods must be a positive integer'))
end

% slice sextupoles
% the number of slices affects the computation of crossing terms.
if nslices == 1
    splitring = ring;
else
    r2=cellfun(@splitelem,ring,'UniformOutput',false);
    splitring=cat(1,r2{:});
end

% prapering the twiss data and magnet data.
% same as computeRDT().
indDQSO=findcells(splitring,'Class','Bend','Quadrupole','Sextupole','Octupole','Multipole');

[~,AVEBETA,AVEMU,AVEDISP,nu,~]=atavedata(splitring,0,1:length(splitring));

s=[0,findspos(splitring,indDQSO)];
betax=AVEBETA(indDQSO,1);
betay=AVEBETA(indDQSO,2);
etax=AVEDISP(indDQSO,1);
phix=AVEMU(indDQSO,1);
phiy=AVEMU(indDQSO,2);
a2L=getcellstruct(splitring,'PolynomA',indDQSO,1,2).*getcellstruct(splitring,'Length',indDQSO);
if any(a2L)
    throw(MException('RDTFluctuation:Unfinished', ...
        'The coupling case cannot be handled in the current version.'))
end
b2L=getcellstruct(splitring,'PolynomB',indDQSO,1,2).*getcellstruct(splitring,'Length',indDQSO);
b2L(isnan(b2L))=0;
b3L=getcellstruct(splitring,'PolynomB',indDQSO,1,3).*getcellstruct(splitring,'Length',indDQSO);
b3L(isnan(b3L))=0;
b4L=getcellstruct(splitring,'PolynomB',indDQSO,1,4).*getcellstruct(splitring,'Length',indDQSO);
b4L(isnan(b4L))=0;
nData=length(indDQSO) + 1;

% calculate the build-up fluctuation of RDTs.
% hjklmqs means h_jklmq(s), build-up RDT fluctuations.
% the accumulated RDT along the longitudinal position.
[h21000s, h30000s, h10110s, h10020s, h10200s, h20001s, h00201s, h10002s,...
          h22000s, h11110s, h00220s,...
          h31000s, h40000s, h20110s, h11200s, h20020s, h20200s, h00310s, h00400s] = RDTbuildupFluct(betax,betay,...
          etax,phix,phiy,b2L,b3L,b4L,nData);

h12000s = conj(h21000s);
h01110s = conj(h10110s);
h01200s = conj(h10020s);
h01020s = conj(h10200s);
          
% one-period RDT
h21000 = h21000s(end);
h30000 = h30000s(end);
h10110 = h10110s(end);
h10020 = h10020s(end);
h10200 = h10200s(end);
h20001 = h20001s(end);
h00201 = h00201s(end);
h10002 = h10002s(end);
h22000 = h22000s(end);
h11110 = h11110s(end);
h00220 = h00220s(end);
h31000 = h31000s(end);
h40000 = h40000s(end);
h20110 = h20110s(end);
h11200 = h11200s(end);
h20020 = h20020s(end);
h20200 = h20200s(end);
h00310 = h00310s(end);
h00400 = h00400s(end);

period_phix = nu(1) * 2 * pi;
period_phiy = nu(2) * 2 * pi;

% calculate natural RDT fluctuation.
% Here fxxxxx is f_xxxxx(0) in (Franchi,2014),
% but fxxxxxs is not equal to f_xxxxx(s); their absolute values are equal.
f21000 = h21000 / (1 - exp(1i * period_phix));
f30000 = h30000 / (1 - exp(1i * 3 * period_phix));
f10110 = h10110 / (1 - exp(1i * period_phix));
f10020 = h10020 / (1 - exp(1i * (period_phix - 2 * period_phiy)));
f10200 = h10200 / (1 - exp(1i * (period_phix + 2 * period_phiy)));
f20001 = h20001 / (1 - exp(1i * 2 * period_phix));
f00201 = h00201 / (1 - exp(1i * 2 * period_phiy));
f10002 = h10002 / (1 - exp(1i * period_phix));

f21000s = f21000 - h21000s;
f30000s = f30000 - h30000s;
f10110s = f10110 - h10110s;
f10020s = f10020 - h10020s;
f10200s = f10200 - h10200s;

% 
f20001 = h20001 / (1 - exp(1i * 2 * period_phix));
f00201 = h00201 / (1 - exp(1i * 2 * period_phiy));
f10002 = h10002 / (1 - exp(1i * period_phix));
f20001s = f20001 - h20001s;
f00201s = f00201 - h00201s;
f10002s = f10002 - h10002s;

f12000 = conj(f21000);
f01110 = conj(f10110);
f01200 = conj(f10020);
f01020 = conj(f10200);
h12000 = conj(h21000);
h01110 = conj(h10110);
h01200 = conj(h10020);
h01020 = conj(h10200);

f31000 = 1i * 6 * (h30000 * f12000 - h12000 * f30000) + h31000;
f40000 = 1i * 3 * (h30000 * f21000 - h21000 * f30000) + h40000;
f20110 = 1i * ((h30000 * f01110 - h01110 * f30000) * 3 ...
               - (h21000 * f10110 - h10110 * f21000) ...
               + (h10200 * f10020 - h10020 * f10200) * 4) + h20110;
f11200 = 1i * ((h10200 * f12000 - h12000 * f10200) * 2 ...
               + (h21000 * f01200 - h01200 * f21000) * 2 ...
               + (h10200 * f01110 - h01110 * f10200) * 2 ...
               + (h10110 * f01200 - h01200 * f10110) * (-2)) + h11200;
f20020 = 1i * (-(h21000 * f10020 - h10020 * f21000) ...
               + (h30000 * f01020 - h01020 * f30000) * 3 ...
               + (h10110 * f10020 - h10020 * f10110) * 2) + h20020;
f20200 = 1i * ((h30000 * f01200 - h01200 * f30000) * 3 ...
               + (h10200 * f21000 - h21000 * f10200) ...
               + (h10110 * f10200 - h10200 * f10110) * (-2)) + h20200;
f00310 = 1i * ((h10200 * f01110 - h01110 * f10200) ...
               + (h10110 * f01200 - h01200 * f10110)) + h00310;
f00400 = 1i * (h10200 * f01200 - h01200 * f10200) + h00400;

f31000 = f31000 / (1 - exp(1i * 2 * period_phix));
f40000 = f40000 / (1 - exp(1i * 4 * period_phix));
f20110 = f20110 / (1 - exp(1i * 2 * period_phix));
f11200 = f11200 / (1 - exp(1i * 2 * period_phiy));
f20020 = f20020 / (1 - exp(1i * (2 * period_phix - 2 * period_phiy)));
f20200 = f20200 / (1 - exp(1i * (2 * period_phix + 2 * period_phiy)));
f00310 = f00310 / (1 - exp(1i * 2 * period_phiy));
f00400 = f00400 / (1 - exp(1i * 4 * period_phiy));

f31000s = (h31000s - (f30000 * h12000s - f12000 * h30000s) * 1i * 6 -f31000);

f40000s = (h40000s - (f30000 * h21000s - f21000 * h30000s) * 1i * 3 -f40000);
f20110s = (h20110s - ((f30000 * h01110s - f01110 * h30000s) * 3 ...
                                - f21000 * h10110s + f10110 * h21000s ...
                                + (f10200 * h10020s - f10020 * h10200s) * 4) * 1i -f20110);
f11200s = (h11200s - (f10200 * (h12000s + h01110s) - f12000 * h10200s ...
                                +f21000 * h01200s - f01200 * (h21000s - h10110s)...
                                -f01110 * h10200s - f10110 * h01200s) * 1i * 2 -f11200);
f20020s = (h20020s - (-f21000 * h10020s + f10020 * (h21000s - h10110s * 2) ...
                                + f30000 * h01020s * 3 - f01020 * h30000s * 3 ...
                                +f10110 * h10020s * 2) * 1i - f20020);
f20200s = (h20200s - (f30000 * h01200s * 3 - f01200 * h30000s * 3 ...
                                + f10200 * (h21000s + h10110s * 2) ...
                                - f21000 * h10200s  ...
                                - f10110 * h10200s * 2) * 1i - f20200);
f00310s = (h00310s - (f10200 * h01110s - f01110 * h10200s ...
                                +f10110 * h01200s - f01200 * h10110s) * 1i - f00310);
f00400s = (h00400s - (f10200 * h01200s - f01200 * h10200s) * 1i - f00400);

natural_fluctuation = struct('s', s, 'f21000',abs(f21000s), 'f30000',abs(f30000s), ...
                    'f10110',abs(f10110s), 'f10020',abs(f10020s), 'f10200',abs(f10200s), ...
                    'f20001',abs(f20001s), 'f00201',abs(f00201s), 'f10002',abs(f10002s), ...
                    'f31000',abs(f31000s), 'f40000',abs(f40000s), 'f20110',abs(f20110s), 'f11200',abs(f11200s), ...
                    'f20020',abs(f20020s), 'f20200',abs(f20200s), 'f00310',abs(f00310s), 'f00400',abs(f00400s));

% calculate build-up fluctuation for multiple periods.
if nperiods == 1
    RDT = struct('h21000', h21000, 'h30000', h30000, 'h10110', h10110, 'h10020', h10020, ...
            'h10200', h10200, 'h20001', h20001, 'h00201', h00201, 'h10002', h10002, ...
            'h22000', h22000, 'h11110', h11110, 'h00220', h00220, ...
            'h31000', h31000, 'h40000', h40000, 'h20110', h20110, 'h11200', h11200, ...
            'h20020', h20020, 'h20200', h20200, 'h00310', h00310, 'h00400', h00400,...
            'dnux_dJx', -4*real(h22000)/pi, 'dnux_dJy', -2*real(h11110)/pi, 'dnuy_dJy',-4*real(h00220)/pi);
    buildup_fluctuation = struct('s', s, 'h21000', h21000s,'h30000', h30000s,...
                    'h10110', h10110s,'h10020', h10020s,'h10200', h10200s,...
                    'h20001', h20001s,'h00201', h00201s,'h10002', h10002s,...
                    'h22000', h22000s, 'h11110', h11110s, 'h00220', h00220s,...
                    'h31000', h31000s,'h40000', h40000s,'h20110', h20110s,'h11200', h11200s,...
                    'h20020', h20020s,'h20200', h20200s,'h00310', h00310s,'h00400', h00400s);
else
    buildup_fluctuation=multiperiod();
    RDT = struct('h21000', buildup_fluctuation.h21000(end), 'h30000', buildup_fluctuation.h30000(end), 'h10110', buildup_fluctuation.h10110(end), 'h10020', buildup_fluctuation.h10020(end), ...
            'h10200', buildup_fluctuation.h10200(end), 'h20001', buildup_fluctuation.h20001(end), 'h00201', buildup_fluctuation.h00201(end), 'h10002', buildup_fluctuation.h10002(end), ...
            'h22000', buildup_fluctuation.h22000(end), 'h11110', buildup_fluctuation.h11110(end), 'h00220', buildup_fluctuation.h00220(end), ...
            'h31000', buildup_fluctuation.h31000(end), 'h40000', buildup_fluctuation.h40000(end), 'h20110', buildup_fluctuation.h20110(end), 'h11200', buildup_fluctuation.h11200(end), ...
            'h20020', buildup_fluctuation.h20020(end), 'h20200', buildup_fluctuation.h20200(end), 'h00310', buildup_fluctuation.h00310(end), 'h00400', buildup_fluctuation.h00400(end),...
            'dnux_dJx', -4*real(buildup_fluctuation.h22000(end))/pi, 'dnux_dJy', -2*real(buildup_fluctuation.h11110(end))/pi, 'dnuy_dJy',-4*real(buildup_fluctuation.h00220(end))/pi);
end

    function newelems=splitelem(elem)
        if isfield(elem,'Length') && elem.Length > 0 ...
                && ~strcmp(elem.PassMethod, 'IdTablePass')...
                && isfield(elem, 'PolynomB') && length(elem.PolynomB) > 2
            newelems=atdivelem(elem,ones(1,nslices)./nslices);
        else
            newelems={elem};
        end
    end

    function multi_period_RDT=multiperiod()
        q21000 = exp(1i * period_phix);
        q30000 = exp(1i * period_phix * 3);
        q10110 = exp(1i * period_phix);
        q10020 = exp(1i * (period_phix - 2 * period_phiy));
        q10200 = exp(1i * (period_phix + 2 * period_phiy));
        q12000 = conj(q21000);
        q01110 = conj(q10110);
        q01200 = conj(q10020);
        q01020 = conj(q10200);
        q03000 = conj(q30000);
        q20001 = exp(1i * 2 * period_phix);
        q00201 = exp(1i * 2 * period_phiy);
        q10002 = exp(1i * period_phix);
        q31000 = exp(1i * 2 * period_phix);
        q40000 = exp(1i * 4 * period_phix);
        q20110 = exp(1i * 2 * period_phix);
        q11200 = exp(1i * 2 * period_phiy);
        q20020 = exp(1i * (2 * period_phix - 2 * period_phiy));
        q20200 = exp(1i * (2 * period_phix + 2 * period_phiy));
        q00310 = exp(1i * 2 * period_phiy);
        q00400 = exp(1i * 4 * period_phiy);

        total_length = nData * nperiods;
        h21000s2 = complex(zeros(total_length, 1));
        h30000s2 = complex(zeros(total_length, 1));
        h10110s2 = complex(zeros(total_length, 1));
        h10020s2 = complex(zeros(total_length, 1));
        h10200s2 = complex(zeros(total_length, 1));
        h20001s2 = complex(zeros(total_length, 1));
        h00201s2 = complex(zeros(total_length, 1));
        h10002s2 = complex(zeros(total_length, 1));
        h31000s2 = complex(zeros(total_length, 1));
        h40000s2 = complex(zeros(total_length, 1));
        h20110s2 = complex(zeros(total_length, 1));
        h11200s2 = complex(zeros(total_length, 1));
        h20020s2 = complex(zeros(total_length, 1));
        h20200s2 = complex(zeros(total_length, 1));
        h00310s2 = complex(zeros(total_length, 1));
        h00400s2 = complex(zeros(total_length, 1));
        s2 = zeros(total_length, 1);
        
        f03000 = conj(f30000);
        f22000s = h22000s + 1i * ((h21000s * f12000 - f21000 * h12000s) * 3 + (h30000s * conj(f30000) - f30000 * conj(h30000s)) * 9);
        f11110s = h11110s + 1i * ((h21000s * f01110 - h01110s * f21000) * 2 - (h12000s * f10110 - h10110s * f12000) * 2 ...
                                - (h10020s * f01200 - h01200s * f10020) * 4 + (h10200s * f01020 - h01020s * f10200) * 4);
        f00220s = h00220s + 1i * ((h10020s * f01200 - h01200s * f10020) ...
                                + (h10200s * f01020 - h01020s * f10200) ...
                                + (h10110s * f01110 - h01110s * f10110));    
        s2(1:nData) = s;
        h21000s2(1:nData) = h21000s;
        h30000s2(1:nData) = h30000s;
        h10110s2(1:nData) = h10110s;
        h10020s2(1:nData) = h10020s;
        h10200s2(1:nData) = h10200s;
        h20001s2(1:nData) = h20001s;
        h00201s2(1:nData) = h00201s;
        h10002s2(1:nData) = h10002s;
        h22000s2(1:nData) = h22000s;
        h11110s2(1:nData) = h11110s;
        h00220s2(1:nData) = h00220s;
        h31000s2(1:nData) = h31000s;
        h40000s2(1:nData) = h40000s;
        h20110s2(1:nData) = h20110s;
        h11200s2(1:nData) = h11200s;
        h20020s2(1:nData) = h20020s;
        h20200s2(1:nData) = h20200s;
        h00310s2(1:nData) = h00310s;
        h00400s2(1:nData) = h00400s;
        for iii = 1:(nperiods-1)
            s2(iii * nData + 1: (iii+1) * nData) = s + findspos(splitring,length(splitring) + 1) * iii;
            h21000s2(iii * nData +1: (iii + 1) * nData) = f21000 + (h21000s - f21000) * q21000^iii;
            h30000s2(iii * nData +1: (iii + 1) * nData) = f30000 + (h30000s - f30000) * q30000^iii;
            h10110s2(iii * nData +1: (iii + 1) * nData) = f10110 + (h10110s - f10110) * q10110^iii;
            h10020s2(iii * nData +1: (iii + 1) * nData) = f10020 + (h10020s - f10020) * q10020^iii;
            h10200s2(iii * nData +1: (iii + 1) * nData) = f10200 + (h10200s - f10200) * q10200^iii;
            h20001s2(iii * nData +1: (iii + 1) * nData) = f20001 + (h20001s - f20001) * q20001^iii;
            h00201s2(iii * nData +1: (iii + 1) * nData) = f00201 + (h00201s - f00201) * q00201^iii;
            h10002s2(iii * nData +1: (iii + 1) * nData) = f10002 + (h10002s - f10002) * q10002^iii;

            h22000s2(iii * nData +1: (iii + 1) * nData) = f22000s(end) * iii +  f22000s ...
                                    +1i * 3 * (f12000 * (f21000 - h21000s) * q21000^iii - f21000 * (f12000 - h12000s) * q12000^iii) ...
                                    +1i * 9 * (f03000 * (f30000 - h30000s) * q30000^iii - f30000 * (f03000 - conj(h30000s)) * q03000^iii);

            h11110s2(iii * nData +1: (iii + 1) * nData) = f11110s(end) * iii + f11110s ...
                                    +1i * 2 * ((f21000 - h21000s) * f01110 * q21000^iii - (f01110 - h01110s) * f21000 * q01110^iii) ...
                                    -1i * 2 * ((f12000 - h12000s) * f10110 * q12000^iii - (f10110 - h10110s) * f12000 * q10110^iii) ...
                                    -1i * 4 * ((f10020 - h10020s) * f01200 * q10020^iii - (f01200 - h01200s) * f10020 * q01200^iii) ...
                                    +1i * 4 * ((f10200 - h10200s) * f01020 * q10200^iii - (f01020 - h01020s) * f10200 * q01020^iii);
            h00220s2(iii * nData +1: (iii + 1) * nData) = f00220s(end) * iii + f00220s ...
                                    +1i * ((f10020 - h10020s) * f01200 * q10020^iii - (f01200 - h01200s) * f10020 * q01200^iii) ...
                                    +1i * ((f10200 - h10200s) * f01020 * q10200^iii - (f01020 - h01020s) * f10200 * q01020^iii) ...
                                    +1i * ((f10110 - h10110s) * f01110 * q10110^iii - (f01110 - h01110s) * f10110 * q01110^iii);   
            h31000s2(iii * nData +1: (iii + 1) * nData) = (f31000 + f31000s * q31000^iii ...
                                    +1i * 6 * f12000 * (f30000 - h30000s) * q30000^iii ...
                                    +1i * 6 * (h12000s - f12000) * f30000 * q12000^iii);
            h40000s2(iii * nData +1: (iii + 1) * nData) = (f40000 + f40000s * q40000^iii ...
                                    + 1i * 3 * f21000 * (f30000 - h30000s) * q30000^iii ...
                                    +1i * 3* (h21000s - f21000) * f30000 * q21000^iii);
            h20110s2(iii * nData +1: (iii + 1) * nData) = (f20110 + f20110s * q20110^iii ...
                                    + (f01110 * (f30000 - h30000s) * 3 * q30000^iii ...
                                    -(f01110 - h01110s) * f30000 * 3 * q01110^iii ...
                                    -f10110 * (f21000 - h21000s) * q21000^iii ...
                                    +(f10110 - h10110s) * f21000 * q10110^iii ...
                                    +f10020 * (f10200 - h10200s) * 4 * q10200^iii ...
                                    - (f10020 - h10020s) * f10200 * 4 * q10020^iii) * 1i);
            h11200s2(iii * nData +1: (iii + 1) * nData) = (f11200 + f11200s * q11200^iii ...
                                    +((f12000 + f01110) * (f10200 - h10200s) * 2 * q10200^iii ...
                                    -(f12000 - h12000s) * f10200 * 2 * q12000^iii ...
                                    +f01200 * (f21000 - h21000s) * 2 * q21000^iii ...
                                    -(f01200 - h01200s) * (f21000 - f10110) * 2 * q01200^iii ...
                                    -(f01110 - h01110s) * f10200 * 2 * q01110^iii ...
                                    +f01200 * (f10110 - h10110s) * (-2) * q10110^iii) * 1i);
            h20020s2(iii * nData +1: (iii + 1) * nData) = (f20020 + f20020s * q20020^iii ...
                                    +(-f10020 * (f21000 - h21000s) * q21000^iii ...
                                    + (f10020 - h10020s) * (f21000 - f10110 * 2) * q10020^iii ...
                                    + f01020 * (f30000 - h30000s) * 3 * q30000^iii ...
                                    - (f01020 - h01020s) * f30000 * 3 * q01020^iii ...
                                    + f10020 * (f10110 - h10110s) * 2 * q10110^iii) * 1i);
            h20200s2(iii * nData +1: (iii + 1) * nData) = (f20200 + f20200s * q20200^iii ...
                                    + (f01200 * (f30000 - h30000s) * 3 * q30000^iii ...
                                    - (f01200 - h01200s) * f30000 * 3 * q01200^iii ...
                                    + (f21000 + 2 * f10110) * (f10200 - h10200s) * q10200^iii ...
                                    - (f21000 - h21000s) * f10200 * q21000^iii ...
                                    + f10200 * (f10110 - h10110s) * (-2) * q10110^iii) * 1i);
            h00310s2(iii * nData +1: (iii + 1) * nData) = (f00310 + f00310s * q00310^iii ...
                                    + (f01110 * (f10200 - h10200s) * q10200^iii ...
                                    - (f01110 - h01110s) * f10200 * q01110^iii ...
                                    + f01200 * (f10110 - h10110s) * q10110^iii ...
                                    - (f01200 - h01200s) * f10110 * q01200^iii) * 1i);
            h00400s2(iii * nData +1: (iii + 1) * nData) = (f00400 + f00400s * q00400^iii ...
                                    + (f01200 * (f10200 - h10200s) * q10200^iii ...
                                    - (f01200 - h01200s) * f10200 * q01200^iii) * 1i);
        end
        multi_period_RDT=struct('s', s2, 'h21000', h21000s2,'h30000', h30000s2,...
                    'h10110', h10110s2,'h10020', h10020s2,'h10200', h10200s2, ...
                    'h20001', h20001s2,'h00201', h00201s2,'h10002', h10002s2,...
                    'h22000', h22000s2,'h11110', h11110s2, 'h00220', h00220s2, ...
                    'h31000', h31000s2,'h40000', h40000s2,'h20110', h20110s2,'h11200', h11200s2,...
                    'h20020', h20020s2,'h20200', h20200s2,'h00310', h00310s2,'h00400', h00400s2);
    end
end

