function bbeat=plot_betabeat(THERING_ref,THERING_mod)
%function plot_betabeat(THERING_ref,THERING_mod)
%
% returns plot of beta beat of THERING_mod respect to THERING_ref

% check that it is the same lattice

indx=1:length(THERING_ref);

if length(THERING_mod) == length(THERING_ref)
    
    Tr=twissring(THERING_ref,0,indx);
    Tm=twissring(THERING_mod,0,indx);
    
    br=cat(1,Tr.beta);
    bm=cat(1,Tm.beta);
    
    bbeat(1)=plot(findspos(THERING_ref,indx),(bm(:,1)-br(:,1))./br(:,1),'.-b');
hold on;
   bbeat(2)=plot(findspos(THERING_ref,indx),(bm(:,2)-br(:,2))./br(:,2),'.-r');

else
    error('Not the same lattice');
end

return


