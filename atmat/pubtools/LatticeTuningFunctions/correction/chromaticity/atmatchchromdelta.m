function [arcchrom0,deltsext]=atmatchchromdelta(arc,c,sxtfams)
% function arcchrom0=atmatchchromdelta(arc,c,sxtfams)
% 
% arc    : at lattice 
% c      : chromaticity to get size(c)=[2,1]
% sxtfams: {[findcells(arc,'FamName','SF1','SF2')],...
%           [findcells(arc,'FamName','SD1','SD2')] }
% 
% adds a common DKsext on the two specified sextupole families 
%
%see also: atmatch

disp('match chromaticity')

variabs=[];

for isextfams=1:length(sxtfams)
  Ksf=cellfun(@(a)a.PolynomB(3),arc(sxtfams{isextfams}),'un',1);
  variabs=[variabs, atVariableBuilder(arc,...
    {@(r,DKsext)setcellstruct(r,'PolynomB',sxtfams{isextfams},Ksf+DKsext,1,3)},...
    {[0]})]; %#ok<*AGROW>
end

% Ksf=cellfun(@(a)a.PolynomB(3),arc(sxtfams{1}),'un',1);
% Ksd=cellfun(@(a)a.PolynomB(3),arc(sxtfams{2}),'un',1);

ConstrChrom=[...
    atlinconstraint(1,{{'chromaticity',{1}}},c(1),c(1),1)...
    atlinconstraint(1,{{'chromaticity',{2}}},c(2),c(2),1)];

tol=1e-4;
[arcchrom0,deltsext]=atmatch(arc,variabs,ConstrChrom,tol,150,0);%,@lsqnonlin);

% tol=1e-8;
% arcchrom0=atmatch(arc,variabs,ConstrChrom,tol,50,3,@lsqnonlin);

% Ksfc=cellfun(@(a)a.PolynomB(3),arcchrom0(sxtfams{1}),'un',1);
% Ksdc=cellfun(@(a)a.PolynomB(3),arcchrom0(sxtfams{2}),'un',1);
% DSF=Ksfc-Ksf;
% DSD=Ksdc-Ksd;

return


