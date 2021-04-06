function result = startsWith(str,pat,varargin)
%startsWith    True if text starts with pattern.
%
%RESULT=startsWith(STR,PAT)
%
%   STR:     Target string: character vector or cell array of character vectors
%            Result is an array with the same size as STR
%   PAT:     Pattern to look for: character vector or cell array of character vectors
%            Result is true is STR starts with any element of PAT
%

if iscell(str)
    result=cellfun(@(s) startsWith(s,pat,varargin{:}), str);
else
    if iscell(pat)
        result = false;
        for p=pat
            if startsWith(str,p{1},varargin{:})
                result=true;
                break;
            endif
        endfor
    else
        ignorecase=getoption(varargin,'IgnoreCase',false);
        patLength = length(pat);
        result=(length(str) >= patLength);
        if result
            if ignorecase
                result=strcmpi(str(1:patLength),pat);
            else
                result=strcmp(str(1:patLength),pat);
            endif
        endif
    endif
endif
endfunction
