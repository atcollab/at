function tune = intfft(X,varargin);
%INTFFT Calculates the tune from interpolated FFT of the trajectory.
% INTFFT(X) X must be a column vector.
%  If X is a matrix - each column is treated as
%  a separate trajectory
% INTFFT(X,GUESS,DELTA) searches for peaks in the FFT spectrum
%  only within the range (X-DELTA ... X+DELTA) 
%  The same values of GUESS and DELTA are used for all columns of X 


[N,L] = size(X);
if L == 0;
    tune = NaN;
    
    return
end
% apply hanning window
%W = diag(sin(pi*(0:N-1)/(N-1)).^2);

%XFFTABS = abs(fft(W*X));
XFFTABS = abs(fft(X));
%Z = zeros(size(XFFTABS));
if nargin==3
    GUESS = varargin{1};
    DELTA = varargin{2};
%     LR = floor(N*(GUESS-DELTA));
%     UR = ceil(N*(GUESS+DELTA));
%     Z(sub2ind(size(XFFTABS),LR,1:length(LR))) = 1;
%     Z(sub2ind(size(XFFTABS),UR+1,1:length(UR)))= -1;
    
    
    searchrange = floor(N*(GUESS-DELTA)):ceil(N*(GUESS+DELTA));
    
    %[psi_k,k] = max(cumsum(Z).*XFFTABS);
    [psi_k,k] = max(XFFTABS(searchrange,:));
    k=k+floor(N*(GUESS-DELTA))-1; 
    
else 
    [psi_k,k] = max(XFFTABS(1:floor(N/2),:));
end

psi_k_plus  = XFFTABS((k+1)+(0:(L-1))*N);
psi_k_minus = XFFTABS((k-1)+(0:(L-1))*N);

G = psi_k_plus>psi_k_minus;

k_r = k+G;
k_l = k-~G;

psi_l = XFFTABS(k_l+(0:(L-1))*N);
psi_r = XFFTABS(k_r+(0:(L-1))*N);

tune = (k_l-1)/N + atan( psi_r*sin(pi/N)./(psi_l + psi_r*cos(pi/N)))/pi;