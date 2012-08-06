function index = findcells(CELLARRAY, field, varargin)
%FINDCELLS performs a search on MATLAB cell arrays of structures
%   
% INDEX = FINDCELLS(CELLARRAY, 'field') 
%   returns indexes of elements that have a field named 'field'   
%
% INDEX = FINDCELLS(CELLARRAY, 'field', VALUE) 
%   returns indexes of elements whose field 'field'
%   is equal to VALUE1, VALUE2, ... or VALUEN. Where VALUE can either be
%   character strings or a number. If its a character string REGULAR
%   expressions can be used.
%
% Example:
%   findcells(THERING,'Length',0, 0.2);  % will match elements of
%                                          lengths 0 and 0.2
%   findcells(THERING,'FamName','SFA','SDA');
%
% See also GETCELLSTRUCT, SETCELLSTRUCT, REGEXPI

index=reshape(find(atgetcells(CELLARRAY,field,varargin{:})),1,[]);
