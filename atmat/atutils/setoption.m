function setoption(name,value)
%SETOPTION	Set AT preference values
%
%SETOPTION('KEY',DEFAULT)
%   Set the default value for the given KEY to DEFAULT. It is an error to set
%   a default for a non-existing KEY. See ATINIPREFS for the list of
%   predefined keys.
%
% KEY:      Key name
% DEFAULT:  New default value for the key
%
%SETOPTION('KEY') Resets the default value for KEY to its inital setting
%
%SETOPTION()      Resets all values to the initial setting
%
%See also GETOPTION, INITOPTIONS

opt=atoptions.instance();
if nargin >= 2
    opt.(name)=value;
elseif nargin >= 1
    opt.reset(name);
else
    opt.instance('reset');
end
