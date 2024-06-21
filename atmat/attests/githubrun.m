% Private. This script runs the AT test suite in a GitHub action. There is
% normally no reason to run it in a user workflow.

tsuite=fullfile(atroot,'attests');
if ~(ispc || ismac)
    % Solution for
    % "Intel MKL ERROR: Parameter 11 was incorrect on entry to DSBEVD."
    % found there:
    % https://fr.mathworks.com/matlabcentral/answers/1983899-intel-mkl-error-when-callying-scipy-from-matlab
    py.sys.setdlopenflags(int32(bitor(int64(py.os.RTLD_LAZY), int64(py.os.RTLD_DEEPBIND))));
end
v=assertSuccess(run(testsuite(tsuite)));
disp(table(v));
