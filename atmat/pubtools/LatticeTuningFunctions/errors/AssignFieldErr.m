function [r,Pnew,Panew]=AssignFieldErr(r,refpos,N,rho,BNn,ANn)
% function r=AssignFieldErr(r,refpos,N,rho,BNn,ANn)
% 
% r : at lattice
% refpos : index of elements for wich the field error has to be applied
% N main component (2=quadrupole,3=sextupole)
% rho: reference radius
%
% the field errors are defined by the magnet designer with the multipole
% expansion
% B=B_N(rho0)*sum_(n=1)^(infty)((BNn+i*ANn)*(z/rho0)^(n-1))
%
% and by AT as
% B=Brho*sum_(n=1)^(infty)((b_n+i*a_n)*(z)^(n-1))
% 
% the input BNn and ANn are the normal and skew field errors at rho0. 
%    as defined by:
% B^(N)_(n) = radius^(n-N)*b_n/b_N
%
% b_n=radius^(N-n)*b_N*B^(N)_(n)
%
% optional output ,Pnew,Panew the PolynomB and PolynomA set in AT.
% 


if nargin<6
    ANn=BNn*0;
end
P=getcellstruct(r,'PolynomB',refpos,1,N);
Pa=getcellstruct(r,'PolynomA',refpos,1,N);

if N==1 && P(1)==0
    % is a dipole.
    P=getcellstruct(r,'BendingAngle',refpos);
end

if length(rho)==1
    Pnew=repmat(rho.^(N-[1:length(BNn)]),length(P),1).*repmat(P,1,length(BNn)).*repmat(BNn,length(P),1);
    Panew=repmat(rho.^(N-[1:length(ANn)]),length(Pa),1).*repmat(P,1,length(ANn)).*repmat(ANn,length(Pa),1); % refer to normal component also skew
elseif length(rho)==size(P,1)
    RNB=[];
    RNA=[];
    
    for irrh=1:size(P,1)
        RNB(irrh,:)=rho(irrh).^(N-[1:length(BNn)]);
        RNA(irrh,:)=rho(irrh).^(N-[1:length(ANn)]);
    end
     
    Pnew=RNB.*repmat(P,1,length(BNn)).*repmat(BNn,length(P),1);
    Panew=RNA.*repmat(P,1,length(ANn)).*repmat(ANn,length(Pa),1);
end

for ir=1:length(refpos)
    pbold=r{refpos(ir)}.PolynomB;
    paold=r{refpos(ir)}.PolynomA;
   
    if ir==1 && length(pbold)>length(Pnew) % pad with zeros
        %disp('padding')
        Pnew=[Pnew zeros(length(refpos),length(pbold)-length(Pnew))];
        Panew=[Panew zeros(length(refpos),length(pbold)-length(Panew))];
    end
    
%     size(pbold)
%     size(paold)
%     size(Pnew)
%     size(Panew)

    Pnew(ir,1:length(pbold))=Pnew(ir,1:length(pbold))+pbold;
    Panew(ir,1:length(paold))=Panew(ir,1:length(paold))+paold;
    
    
    r{refpos(ir)}.PolynomB=Pnew(ir,:);
    r{refpos(ir)}.PolynomA=Panew(ir,:); % update both A and B fields! 
    r{refpos(ir)}.MaxOrder=size(Pnew,2)-1; % -1 for c-integrators consitency
end


return