function b=bmat(ring)
%construct the B matrix that transforms into variables including dispersion.

[ld,prm]=atx(ring);

%disp = cat(2,ld.dispersion);
disp = ld(1).Dispersion;

%b = [1 0 0 0 0 -disp(1); 0 1 0 0 0 -disp(2); 0 0 1 0 0 -disp(3); 0 0 0 1 0 -disp(4); disp(2) -disp(1) disp(4) -disp(3) 1 0; 0 0 0 0 0 1];
b = [1 0 0 0 -disp(1) 0; 0 1 0 0 -disp(2) 0; 0 0 1 0 -disp(3) 0; 0 0 0 1 -disp(4) 0; 0 0 0 0 1 0; disp(2) -disp(1) disp(4) -disp(3) 0 1];