function hjklm=calc_hjklm_Wang(ring,j,k,l,m)
%dnudx=nuamp_analyt(ring,sym)

if nargin==1
   sym=1;
end

%here we implement the formulas of Wang ChunXi for higher order
%driving terms

%[td, tune,chrom] = twissring(ring,0,1:length(ring), 'chrom', 1e-5);
[lindat,nu,xi]=atlinopt(ring,0,1:length(ring));

betaxy = cat(1, lindat.beta);
betax = betaxy(:,1);
betay = betaxy(:,2);

mu = cat(1, lindat.mu);
mux = mu(:,1);
muy = mu(:,2);

 hjklm=0;

 
%now we need to add up h_jklm contribution from the sextupoles.
% 
%basically, we define
%hjklm = sum_u b3L_u beta_x^(j+k)/2 beta_y^(l+m)/2
%e^(i(j-k)phi_x+i(l-m)phi_y)
 
index=findcells(ring,'PolynomB');

for q=1:length(ring)
  for qq=1:length(ring)
    b3l_q=ring{index(q)}.PolynomB(3);
    b3l_qq=ring{index(qq)}.PolynomB(3);
    hjklm=hjklm+b3l_q*b3l_qq*
    ...
    betax(index(q))^((j+k)/2)*betay(index(q))^((j+k)/2)*exp(i*((j-k)*mux(index(q))+(l-m)*muy(index(q))));
end



