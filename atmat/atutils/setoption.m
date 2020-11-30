function setoption(varargin)
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
%SETOPTION() Resets all predefined keys to their initial setting.
%
%See also GETOPTION, INITOPTIONS

atprefutil('set',varargin{:});
end
