function [t, errorstr] = isatelem(ELEM,varargin)
%ISATELEM tests if an input argument is a valid AT element.
% 
%  A valid AT element is a MATLAB structure with required 
%   fields 'Length', 'PassMethod', and a set of data fields,
%   specific to the PassMethod used.
%   
%  [TEST, ERRORSTR] = ISATELEM(ELEM)
%                   = ISATELEM(ELEM, 'display')
%
%  TEST     - test result,  1 = valid AT element
%  ERRORSTR - multi-line error message
%    
%  See also PASSMETHOD, ATELEM


errorstr = [];

if ~isstruct(ELEM)
    errorstr = [errorstr,sprintf('%s\n','Input is not a MATLAB structure')];
else
    if ~isfield(ELEM,'PassMethod');
        errorstr = [errorstr,sprintf('%s\n','Required field ''PassMethod'' is missing')];
    else % check if ELEM has all fields required by PassMethod function
        EXISTRESULT = exist(ELEM.PassMethod);
        if EXISTRESULT == 3
            
            
            try % Try to propagate a test particle
                temp = feval(ELEM.PassMethod,ELEM, [0 0 0 0 0 0]');
                
            catch
                errorstr = [errorstr,sprintf('%s\n',['Specified PassMethod m-file: ''',...
                        (ELEM.PassMethod), ''' returned an error'])];
            end
            
            ReqFields = feval(ELEM.PassMethod);
            
            for field = 1:length(ReqFields)
                if ~isfield(ELEM,ReqFields{field})
                    errorstr = [errorstr,sprintf('%s\n',['Required field ''',ReqFields{field}...
                                ,''' is missing'])];
                end
            end
            
            
            
        elseif EXISTRESULT == 2
            
                       
            try % Try to propagate a test particle
                temp = feval(ELEM.PassMethod,ELEM, [0 0 0 0 0 0]');
                
            catch
                errorstr = [errorstr,sprintf('%s\n',['Specified PassMethod m-file: ''',...
                        (ELEM.PassMethod), ''' returned an error'])];
            end
            
        else
            errorstr = [errorstr,sprintf('%s\n',['Specified PassMethod mex-file or m-file: ''',...
                        (ELEM.PassMethod), '.',mexext,''' does not exist'])];
        end
    end
    
end


if isempty(errorstr)
    t = 1;
else
    t = 0;
end

if any(strncmpi(varargin,'disp',4))
    disp(errorstr);
end