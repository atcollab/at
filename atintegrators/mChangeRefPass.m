function [Out1, Out2]  = mChangeRefPass(Elem,Rin)
%MCHANGEREFPASS Change the reference energy by scaling

if nargin
    scaling=Elem.FieldScaling;
    Out1 = Rin;
    Out1(2,:)=Rin(2,:)/scaling;
    Out1(4,:)=Rin(4,:)/scaling;
    Out1(5,:)=(Rin(5,:) + (1.0-scaling))/scaling;
%     Out1(5,:)=(Rin(5,:) + (1.0-scaling));
else % If called with no input args - return lists of required and optional fields
    Out1 = {'FieldScaling'};
    Out2 = {};
end
