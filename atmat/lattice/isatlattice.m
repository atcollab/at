function [T, badlist, errorstr] = isatlattice(LATTICE,varargin)
%ISATLATTICE tests if an input argument is a valid AT lattice.
% 
%  A valid AT lattice is a MATLAB cell array of valid AT elements
%   
%  [TEST, BADLIST, ERRORSTR] = ISATLATTICE(ELEM, ['option1',...])
%
%  Allowed otions: 
%   'display' - print problems to command window;
%   'stop'    - return after the first problem is found
% 
%  TEST     - test result,  1 = valid AT element
%  ERRORSTR - multi-line error message
%    
%  See also ISATELEM, ATELEM, ATLATTICE

T = 1;
errorstr = [];
badlist = [];

DISPLAYFLAG = 0;
STOPFLAG = 0;

if any(strncmpi(varargin,'disp',4))
    DISPLAYFLAG = 1;
end

if any(strncmpi(varargin,'stop',4))
    STOPFLAG = 1;
end

if ~ iscell(LATTICE)
    errorstr = [errorstr,sprintf('%s\n','Input is not a MATLAB cell array')];
    T = 0;
else
    for k = 1:length(LATTICE)
        [t, s] = isatelem(LATTICE{k});
        if ~t
            badlist = [badlist, k];
            T = 0;
            if DISPLAYFLAG
                s = [sprintf('\nInvalid lattice element # %s\n', int2str(k)),s];
                disp(s);
            end
            if STOPFLAG
                return
            end
        end
    end
    errorstr = 'One or more elements are not valid AT elements';
end
