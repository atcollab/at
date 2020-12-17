function [nu,chi]=atdampingrates(m66)
%ATDAMPINGRATES find tunes and damping rates from one map matrix with radiation
%
%[NU,CHI]=ATDAMPINGRATES(M66)
%
%note that in order to find the damping times, one needs the revolution
%time T0, then
%tau1 = T0/chi1, tau2 = T0/chi2, tau3 = T0/chi3

aa=amat(m66);

Rmat=aa\m66*aa;

[chi,nu]=cellfun(@decode,num2cell(reshape(1:size(m66,1),2,[]),1));

    function [chi,nu]=decode(range)
        matr=Rmat(range,range);
        % nu=mod(atan2(matr(1,2)-matr(2,1),matr(1,1)+matr(2,2))/2/pi,1);
        % chi=-log(sqrt(det(matr)));
        ev=eig(matr);
        ev1log=log(ev(1));
        chi=-real(ev1log);
        nu=mod(imag(ev1log)/2/pi,1);
    end
end
