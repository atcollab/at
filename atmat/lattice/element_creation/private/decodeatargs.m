function [rsrc,varargout] = decodeatargs(default_values,args)
%EXTRACTRES separates arguments and resources
%[RSRC,ARGS]=EXTRACTRES(DEFARGS,ARGLIST)
%   DEFARGS must have length >= 2

nopass=@(arg) ischar(arg) && ~isempty(arg) && isempty(regexp(arg,'.*Pass$','once'));
chararg=find(cellfun(nopass,[args {'x'}]),1);
rsrc=args(chararg:end);
varargout=parseargs(default_values,args(1:chararg-1));
end
