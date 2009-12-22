function [Out1, Out2]  = mdriftpass(DriftData,Rin);
% MDRIFTPASS - example of pass method in matlab
%  Same physics and calling syntax as DriftPass.c

if nargin
    Out1 = Rin;
    NormL = DriftData.Length./(1+Rin(5,:));
    Out1(1,:) = Out1(1,:) + Rin(2,:).*NormL;
    Out1(3,:) = Out1(3,:) + Rin(4,:).*NormL;
    Out1(6,:) = Out1(6,:) + NormL.*(Rin(2,:).*Rin(2,:) + Rin(4,:).*Rin(4,:))./(1+Rin(5,:))/2;
else % If called with no input args - return lists of required and optional fields
    Out1 = {'Length'};
    Out2 = {};

end
