function rrand=seterrorrand(...
    r,...             % nominal lattice
    positions,...      % positions where to apply the error
    errorsetfunct,... % function handle to set the error
    seed,...       % seed [m]
    sigma,...      % sigma [m]
    nsigma,...     % truncation [m]
    exixstingerrval)
%
%  rwave=seterrorrand(...
%     r,...             % nominal lattice
%     position,...      % positions where to apply the error
%     errorsetfunct,... % function handle to set the error
%     seed,...       % seed [m]
%     sigma,...      % sigma [m]
%     nsigma)        % truncation [m]% >0 or set to 2
%
% sets error random. errors are defined and applid by the function
% errorsetfunct with signature r=errorsetfunct(r,positions,erroval)
%
% if seed==0 the random stream is not updated.
%
%
%see also: TruncatedGaussian

if sigma==0
    sigma=1e-15;
end

if nsigma<0
    nsigma=2;
end

if seed~=0
    disp(['Setting Random Stream to seed: ' num2str(seed)]);
    % set seed
    s = RandStream('mcg16807','Seed',seed);
    RandStream.setGlobalStream(s);
else
   % disp('Using previously set random stream')
end

% define vector of errors
errorvalues=TruncatedGaussian(...
    sigma,...
    nsigma*sigma, ...
    size(positions));

if ~isrow(errorvalues)
errorvalues=errorvalues';
end

% default existing error values
if nargin<7
    exixstingerrval=zeros(size(errorvalues));
end

if ~isrow(exixstingerrval)
exixstingerrval=exixstingerrval';
end

% apply error
rrand=errorsetfunct(r,positions,exixstingerrval+errorvalues);

return


