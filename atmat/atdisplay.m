function atdisplay(verbosity,message)
%ATDISPLAY checks the verbosity level in the global variable GLOBVAL
%          and displays message if this is greater than the verbosity
%          for this message.
%
%  See also NUMDIFPARAMS

global GLOBVAL

if isfield(GLOBVAL,'verbosity')&& verbosity<= GLOBVAL.verbosity
    disp(message)
end
