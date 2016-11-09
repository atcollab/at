function rerr=setFieldIntegralError(r0,rerr,indx,order,Nsigma,sigmaperc)
%function rerr=setFieldIntegralError(r0,rerr,indx,order,Nsigma,sigmaperc)
%
%
%see also: ApplyErrorRand
rerr=PadPolynomAB(rerr);
r0=PadPolynomAB(r0);

% original value
try
    kl=atgetfieldvalues(rerr,indx,'PolynomB',{1,order});
    kl(isnan(kl))=0;
catch exc
    kl=zeros(size(dipindx));
end

% design value
if order==1
    kl_0=atgetfieldvalues(r0,indx,'BendingAngle',{1,1});%bending angle
else
    kl_0=atgetfieldvalues(r0,indx,'PolynomB',{1,order});%original model set of polynomB
end

% errors set
kl_err=TruncatedGaussian(sigmaperc,Nsigma*sigmaperc,length(indx))'.*kl_0(:);% compute error starting from model values

% previous set
try
    kle=atgetfieldvalues(rerr,indx,'PolynomBErr',{1,order});
    kle(isnan(kle))=0;
catch exc
    kle=zeros(size(indx));
end

newkl=kl(:)-kle(:)+kl_err(:);

% new set
rerr=atsetfieldvalues(rerr,indx,'PolynomB',{1,order},newkl);
rerr=atsetfieldvalues(rerr,indx,'PolynomBErr',{1,order},kl_err);
rerr=atsetfieldvalues(rerr,indx,'PolynomB0',{1,order},kl);

return