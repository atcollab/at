function varargout=at_gpuinfo() %#ok<STOUT>
% INFO = at_gpuinfo()
% INFO          1xn structure with the following fields:
%               Name: GPU name
%               Version: CUDA compute capability (? for OpenCL)
%               CoreNumber: Multi processor number
%               Platform: Platform name
error('at:missingMex','missing MEX file.');
end
