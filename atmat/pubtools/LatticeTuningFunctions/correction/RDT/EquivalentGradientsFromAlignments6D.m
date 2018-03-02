function [kn,ks,ind]=EquivalentGradientsFromAlignments6D(r,inCOD)
%EQUIVALENTGRADIENTSFROMALIGNMENTS6D Estimated normal quad gradients from sext offsets
%[kn,    1) estimated normal quad gradients from sext offsets, quad
%           errors in quadrupoles and sextupoles.
% ks,    2) estimated skew quad gradients from sext offsets, quad
%           errors in quadrupoles, quadrupole rotation.
% ind    3) indexes of locations at wich kn and ks are found
% ]=EquivalentGradientsFromAlignments6D(
% r,     1) AT lattice structure with errors
% inCOD
% )
%
% the function finds the closed orbit at sextupoles and converts it to
% equivalent quadrupole and skew quadrupole gradients for the computation
% of skew and normal quadrupole RDT 
% quadrupole rotations are also converted in skew quadrupole gradients.
% 
% it returns the complete list of normal (kn) and skew (ks) quadrupole
% gradients at the given indexes (ind) ( not integrated, PolynomB)
% 
% 
%  See also

% quadrupole and skew quadrupole errors are introduced via COD in
% sextupoles
indsext=find(atgetcells(r,'Class','Sextupole'))';

b3=atgetfieldvalues(r,indsext,'PolynomB',{1,3});

oin=findorbit6(r,indsext,inCOD);% orbit at entrance of sextupole DO NOT USE HERE findorbit6Err! 
oout=findorbit6(r,indsext+1,inCOD); % orbit at exit of sextupole

xmisal=cellfun(@(a)getT1(a,1),r(indsext));
ymisal=cellfun(@(a)getT1(a,3),r(indsext));

Dx=(oout(1,:)+oin(1,:))/2; % orbit average in sextupole
Dy=(oout(3,:)+oin(3,:))/2; % orbit average in sextupole

% quarupole errors in sextupoles
kn_sext_err=cellfun(@(a)a.PolynomB(2),r(indsext));
ks_sext_err=cellfun(@(a)a.PolynomA(2),r(indsext));

kn_sext=-2.*b3.*(-Dx+xmisal')'+kn_sext_err;
ks_sext=-2.*b3.*(-Dy+ymisal')'+ks_sext_err;

% quadrupole rotations
indquad=find(atgetcells(r,'Class','Quadrupole'))';

kn2=atgetfieldvalues(r,indquad,'PolynomB',{1,2});
ks2=atgetfieldvalues(r,indquad,'PolynomA',{1,2});
srot=cellfun(@(a)getR1(a),r(indquad));

kn_quad=(1-srot).*kn2;
ks_quad=-srot.*kn2+ks2;

% all elements with PolynomB, not sextupoles or quadrupoles
indPolB=find(atgetcells(r,'PolynomB') & ...
    ~atgetcells(r,'Class','Quadrupole') & ...
    ~atgetcells(r,'Class','Sextupole' ))';

NpolB=cellfun(@(a)length(a.PolynomB),r(indPolB));
NpolA=cellfun(@(a)length(a.PolynomA),r(indPolB));

indPolB=indPolB(NpolB>=2 & NpolA>=2 & ~ismember(indPolB,[indquad indsext])');

kn_all=cellfun(@(a)a.PolynomB(2),r(indPolB));%-kn0_all;
ks_all=cellfun(@(a)a.PolynomA(2),r(indPolB));%-ks0_all;

[ind,ord]=sort([indsext,indquad,indPolB]);

kn=[kn_sext;kn_quad;kn_all];
ks=[ks_sext;ks_quad;ks_all];

% integrated strengths
%L=cellfun(@(a)a.Length,r(ind));

kn=kn(ord);%.*L;
ks=ks(ord);%.*L;


return


function t1=getT1(a,ind)
t1=0;
if isfield(a,'T1')
    t1=-a.T1(ind);
end
return


function r1=getR1(a)
r1=0;
if isfield(a,'R1')
    r1=asin(a.R1(1,3));
end
return
