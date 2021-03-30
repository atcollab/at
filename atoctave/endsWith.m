function result = endsWith(str,pat,varargin)
%endsWith    True if text ends with pattern.
%
%RESULT=endsWith(STR,PAT)
%
%   STR:     Target string: character vector or cell array of character vectors
%            Result is an array with the same size as STR
%   PAT:     Pattern to look for: character vector or cell array of character vectors
%            Result is true is STR ends with any element of PAT
%

if iscell(str)
    result=cellfun(@(s) endsWith(s,pat,varargin{:}), str);
else
    if iscell(pat)
        result = false;
        for p=pat
            if endsWith(str,p{1},varargin{:})
                result=true;
                break;
            end
        endfor
    else
        ignorecase=getoption(varargin,'IgnoreCase',false);
        patLength = length(pat);
        result=(length(str) >= patLength);
        if result
            if ignorecase
                result=strcmpi(str(end-patLength+1:end),pat);
            else
                result=strcmp(str(end-patLength+1:end),pat);
            endif
        endif
    endif
endif
endfunction
