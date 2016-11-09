function orbit = findorbit6Err(RING, indbpm,varargin)
% findorbit6 with bpm reading errors
%
%see also findorbit6 bpm_process bpm_matrices

orbit = findorbit6(RING, indbpm, varargin{:});

% BPM errors used only if class is Monitor
%useind=ismember(indbpm,find(atgetcells(ring,'Class','Monitor'))); 
useind=1:length(indbpm);

[rel,tel,trand] = bpm_matrices(RING(indbpm(useind)));
bpmreading = bpm_process(orbit([1,3],:),rel,tel,trand);
orbit(1,:)=bpmreading(1,:);
orbit(3,:)=bpmreading(2,:);

return
