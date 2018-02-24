function orbit = findorbit4Err(RING, dP, indbpm, varargin)
%findorbit4Err - 4x4 closed orbit with BPM errors
%
%  INPUTS
%    1. RING Ring structure
%    2. dP   Energy offset
%    3. indbpm Indices of BPMs (Extended BPM any point to monitor)
%
%  OUTPUTS
%    1. orbit 4x4 orbit in meters
%  
%  NOTES 
%   1. Alignment errors (T1,T2) are considered here and the BPM reading is
%      modified accordingly.
% 
%  See also findorbit4 ApplyBPMErr bpm_matrices

% Get 4x4 closed orbit
orbit = findorbit4(RING, dP,indbpm, varargin{:});

% Get transformation matrices for applying error on BPM readings
[rel,tel,trand] = bpm_matrices(RING(indbpm));

% Process BPM reading according to alignment errors
bpmreading = bpm_process(orbit([1,3],:),rel,tel,trand);

orbit(1,:)=bpmreading(1,:);
orbit(3,:)=bpmreading(2,:);


return
