% Private. This script runs the AT test suite in a GitHub action. There is
% normally no reason to run it in a user workflow.

tsuite=fullfile(atroot,'attests');
% if ispc || ismac
%     v=assertSuccess(run(testsuite(tsuite)));
% else
%     v=assertSuccess(run(testsuite(tsuite,'Tag','GitHub')));
% end
if ~(ispc || ismac)
    py.sys.setdlopenflags(int32(bitor(int64(py.os.RTLD_LAZY), int64(py.os.RTLD_DEEPBIND))));
end
v=assertSuccess(run(testsuite(tsuite)));
disp(table(v));
