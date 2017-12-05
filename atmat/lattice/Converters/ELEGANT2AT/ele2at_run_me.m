%test_elegant_converter

elegant2at('test.new',3);

load test_AT.mat;

atplot(THERING);

try
    plotbeta;
catch
    disp('elegant 2 at conversion unsuccessful. Check if files are in the right folder');
end
