function orbit = findorbit4Err(RING, D, indbpm,varargin)
% findorbit4 with bpm errors
% 
% also alignment errors (T1,T2) are considered here and the BPM reading is
% modified accordingly.
% 
%see also findorbit4 ApplyBPMErr

orbit = findorbit4(RING, D,indbpm, varargin{:});

[rel,tel,trand] = bpm_matrices(RING(indbpm));
bpmreading = bpm_process(orbit([1,3],:),rel,tel,trand);
orbit(1,:)=bpmreading(1,:);
orbit(3,:)=bpmreading(2,:);


return
