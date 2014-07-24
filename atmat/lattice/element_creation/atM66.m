function elem=atM66(fname,m66)
%   atm66 creates an element that applies an arbitrary matrix m66

elem=atbaselem(fname,'Matrix66Pass','M66',m66,'Class','Matrix66','Length',0.0);