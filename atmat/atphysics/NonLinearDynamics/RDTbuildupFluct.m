function [h21000s, h30000s, h10110s, h10020s, h10200s, h20001s, h00201s, h10002s,...
          h22000s, h11110s, h00220s,...
          h31000s, h40000s, h20110s, h11200s, h20020s, h20200s, h00310s, h00400s] = RDTbuildupFluct(betax,betay,...
          etax,phix,phiy,b2L,b3L,b4L,nData)
%BUILDUPRDTFLUCTUATION build-up fluctuation of RDTs
%
%   This function can show the build-up and cancellation of RDTs.
%   DON'T use this function directly!!! 
%   use computeRDTfluctuation()

h21000s = complex(zeros(nData, 1));
h30000s = complex(zeros(nData, 1));
h10110s = complex(zeros(nData, 1));
h10020s = complex(zeros(nData, 1));
h10200s = complex(zeros(nData, 1));
h20001s = complex(zeros(nData, 1));
h00201s = complex(zeros(nData, 1));
h10002s = complex(zeros(nData, 1));
h22000s = complex(zeros(nData, 1));
h11110s = complex(zeros(nData, 1));
h00220s = complex(zeros(nData, 1));
h31000s = complex(zeros(nData, 1));
h40000s = complex(zeros(nData, 1));
h20110s = complex(zeros(nData, 1));
h11200s = complex(zeros(nData, 1));
h20020s = complex(zeros(nData, 1));
h20200s = complex(zeros(nData, 1));
h00310s = complex(zeros(nData, 1));
h00400s = complex(zeros(nData, 1));

h21000 = 0;
h30000 = 0;
h10110 = 0;
h10020 = 0;
h10200 = 0;
h20001 = 0;
h00201 = 0;
h10002 = 0;
h22000 = 0;
h11110 = 0;
h00220 = 0;
h31000 = 0;
h40000 = 0;
h20110 = 0;
h11200 = 0;
h20020 = 0;
h20200 = 0;
h00310 = 0;
h00400 = 0;
for ii=1:(nData-1)
    h21000s(ii) = h21000;
    h30000s(ii) = h30000;
    h10110s(ii) = h10110;
    h10020s(ii) = h10020;
    h10200s(ii) = h10200;
    h20001s(ii) = h20001;
    h00201s(ii) = h00201;
    h10002s(ii) = h10002;
    h22000s(ii) = h22000;
    h11110s(ii) = h11110;
    h00220s(ii) = h00220;
    h31000s(ii) = h31000;
    h40000s(ii) = h40000;
    h20110s(ii) = h20110;
    h11200s(ii) = h11200;
    h20020s(ii) = h20020;
    h20200s(ii) = h20200;
    h00310s(ii) = h00310;
    h00400s(ii) = h00400;

    betax_i = betax(ii);
    betay_i = betay(ii);
    etax_i = etax(ii);
    phix_i = phix(ii);
    phiy_i = phiy(ii);
    if b2L(ii)~=0
        b2l = b2L(ii);
        h20001 = h20001 + b2l * betax_i * exp(0 + 2i * phix_i) / 8;
        h00201 = h00201 - b2l * betay_i * exp(0 + 2i * phiy_i) / 8;
        h10002 = h10002 + b2l * sqrt(betax_i) * etax_i * exp(0 + 1i * phix_i) / 2;
    end
    if b3L(ii)~=0
        b3l = b3L(ii);
        h20001 = h20001 - b3l * betax_i * etax_i * exp(0 + 1i * 2 * phix_i) / 4;
        h00201 = h00201 + b3l * betay_i * etax_i * exp(0 + 1i * 2 * phiy_i) / 4;
        h10002 = h10002 - b3l * sqrt(betax_i) * etax_i^2 * exp(0 + 1i * phix_i) / 2;

        h21000j = -b3l * betax_i^1.5 * exp(1i * phix_i) / 8;
        h30000j = -b3l * betax_i^1.5 * exp(1i * 3 * phix_i) / 24;
        h10110j = b3l * sqrt(betax_i) * betay_i * exp(1i * phix_i) / 4;
        h10020j = b3l * sqrt(betax_i) * betay_i * exp(1i * (phix_i - 2 * phiy_i)) / 8;
        h10200j = b3l * sqrt(betax_i) * betay_i * exp(1i * (phix_i + 2 * phiy_i)) / 8;

        h12000j = conj(h21000j);
        h01110j = conj(h10110j);
        h01200j = conj(h10020j);

        h12000 = conj(h21000);
        h01110 = conj(h10110);
        h01200 = conj(h10020);

        h22000 = h22000 + 1i * ((h21000 * h12000j - h12000 * h21000j) * 3 ...
                  + (h30000 * conj(h30000j) - conj(h30000) * h30000j) * 9);

        h11110 = h11110 + 1i * ((h21000 * h01110j - h01110 * h21000j) * 2 ...
                  - (h12000 * h10110j - h10110 * h12000j) * 2 ...
                  - (h10020 * h01200j - h01200 * h10020j) * 4 ...
                  + (h10200 * conj(h10200j) - conj(h10200) * h10200j) * 4);

        h00220 = h00220 + 1i * ((h10020 * h01200j - h01200 * h10020j) ...
                  + (h10200 * conj(h10200j) - conj(h10200) * h10200j) ...
                  + (h10110 * h01110j - h01110 * h10110j));

        h31000 = h31000 + 1i * (h30000 * h12000j - h12000 * h30000j) * 6;

        h40000 = h40000 + 1i * (h30000 * h21000j - h21000 * h30000j) * 3;

        h20110 = h20110 + 1i * ((h30000 * h01110j - h01110 * h30000j) * 3 ...
                  - (h21000 * h10110j - h10110 * h21000j) ...
                  + (h10200 * h10020j - h10020 * h10200j) * 4);

        h11200 = h11200 + 1i * ((h10200 * h12000j - h12000 * h10200j) * 2 ...
                  + (h21000 * h01200j - h01200 * h21000j) * 2 ...
                  + (h10200 * h01110j - h01110 * h10200j) * 2 ...
                  - (h10110 * h01200j - h01200 * h10110j) * 2);

        h20020 = h20020 + 1i * (-(h21000 * h10020j - h10020 * h21000j) ...
                  + (h30000 * conj(h10200j) - conj(h10200) * h30000j) * 3 ...
                  + (h10110 * h10020j - h10020 * h10110j) * 2);

        h20200 = h20200 + 1i * ((h30000 * h01200j - h01200 * h30000j) * 3 ...
                  + (h10200 * h21000j - h21000 * h10200j) ...
                  - (h10110 * h10200j - h10200 * h10110j) * 2);

        h00310 = h00310 + 1i * ((h10200 * h01110j - h01110 * h10200j) ...
                  + (h10110 * h01200j - h01200 * h10110j));

        h00400 = h00400 + 1i * (h10200 * h01200j - h01200 * h10200j);

        h21000 = h21000 + h21000j;
        h30000 = h30000 + h30000j;
        h10110 = h10110 + h10110j;
        h10020 = h10020 + h10020j;
        h10200 = h10200 + h10200j;
    end
    if b4L(ii)~=0
        b4l = b4L(ii);
        h22000 = h22000 - 3 * b4l * betax_i^2 / 32;
        h11110 = h11110 + 3 * b4l * betax_i * betay_i / 8;
        h00220 = h00220 - 3 * b4l * betay_i^2 / 32;
        
        h31000 = h31000 - b4l * betax_i^2 * exp(1i * 2 * phix_i) / 16;
        h40000 = h40000 - b4l * betax_i^2 * exp(1i * 4 * phix_i) / 64;
        h20110 = h20110 + 3 * b4l * betax_i * betay_i * exp(1i * 2 * phix_i) / 16;
        h11200 = h11200 + 3 * b4l * betax_i * betay_i * exp(1i * 2 * phiy_i) / 16;
        h20020 = h20020 + 3 * b4l * betax_i * betay_i * exp(1i * (2 * phix_i - 2 * phiy_i)) / 32;
        h20200 = h20200 + 3 * b4l * betax_i * betay_i * exp(1i * (2 * phix_i + 2 * phiy_i)) / 32;
        h00310 = h00310 - b4l * betay_i^2 * exp(1i * 2 * phiy_i) / 16;
        h00400 = h00400 - b4l * betay_i^2 * exp(1i * 4 * phiy_i) / 64;
    end
end
    h21000s(nData) = h21000;
    h30000s(nData) = h30000;
    h10110s(nData) = h10110;
    h10020s(nData) = h10020;
    h10200s(nData) = h10200;
    h20001s(nData) = h20001;
    h00201s(nData) = h00201;
    h10002s(nData) = h10002;
    h22000s(nData) = h22000;
    h11110s(nData) = h11110;
    h00220s(nData) = h00220;
    h31000s(nData) = h31000;
    h40000s(nData) = h40000;
    h20110s(nData) = h20110;
    h11200s(nData) = h11200;
    h20020s(nData) = h20020;
    h20200s(nData) = h20200;
    h00310s(nData) = h00310;
    h00400s(nData) = h00400;
end

