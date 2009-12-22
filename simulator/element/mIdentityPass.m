function [Out1, Out2] = midentitypass(IdentityStruct, Rin);
% MDRIFTPASS - example of pass method in matlab
%  the samephysics as DriftPass.c

if nargin
    Out1 = Rin;
else % If called with no input args - return lists of required and optional fields
    Out1 = {};
    Out2 = {};
end


