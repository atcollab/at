function naff_example
% Example to test naff within matlab

NT = 9996; % divided by 6

% simple quasiperiodic (even period) motion 
y =2+0.1*cos(pi*(0:NT-1))+0.00125*cos(pi/3*(0:NT-1));
yp=2+0.1*sin(pi*(0:NT-1))+0.00125*sin(pi/3*(0:NT-1));

[nu ampl phase] = calcnaff(y,yp,1,'Debug'); % with windowing

str = [
'NFS = 3\n' ...
'AMPL= 2.000000e+00+i* 2.000000e+00 abs(AMPL)= 2.828427e+00 arg(AMPL)= 7.853982e-01 FREQ= 3.904977e-17\n'...
'AMPL= 1.000000e-01+i* 5.638068e-13 abs(AMPL)= 1.000000e-01 arg(AMPL)= 5.638068e-12 FREQ=-3.141593e+00\n' ...
'AMPL= 1.250000e-03+i*-4.767346e-14 abs(AMPL)= 1.250000e-03 arg(AMPL)=-3.813877e-11 FREQ= 1.047198e+00\n' ...
];
fprintf('* one should get the following *\n');
fprintf(str);