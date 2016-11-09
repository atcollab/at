function rwave=seterrorwave(...
    r,...             % nominal lattice
    positions,...      % positions where to apply the error
    errorsetfunct,... % function handle to set the error
    wavelength,...    % wavelength [m]
    amplitude,...     % amplitude [m]
    existingerr)        
%
%  rwave=seterrorwave(...
%     r,...             % nominal lattice
%     position,...      % positions where to apply the error
%     errorsetfunct,... % function handle to set the error
%     wavelength,...    % array of wavelengths [m]
%     amplitude)        % array of amplitudes [m]
%
% sets error waves. errors are defined and applid by the function
% errorsetfunct with signature r=errorsetfunct(r,positions,erroval)
% 
%see also: findspos
if nargin<6
existingerr=zeros(size(positions));
end

nw=length(amplitude);
np=length(positions);

spos=findspos(r,positions);

%errorvalues=zeros(size(positions));

sposdivwave = bsxfun(@rdivide, repmat(spos,nw,1), wavelength');
sinsposdivwave=sin(2*pi*sposdivwave);

sinusoidsmat=repmat(amplitude,np,1).*sinsposdivwave';
errorvalues=sum(sinusoidsmat,2)';

%plot(spos,errorvalues)
%errorvalues=amplitude*sin(2*pi*spos/wavelength);

rwave=errorsetfunct(r,positions,existingerr+errorvalues);

return