function arctune0=atmatchtunedelta(arc,tune,quadfams)
% function arcchrom0=atmatchtunedelta(arc,c,quadfams)
% 
% arc    : at lattice 
% tune      : tune to get (with integer part) size(tune)=[2,1]
% quadfams: {[findcells(arc,'FamName','QF1','QF2')],...
%          [findcells(arc,'FamName','QD1','QD2')] }
% 
% delta on quadrupole families 
%
% fits the tune to the desired values, including the integer part.

%disp('match tunes')

variabs=[];

for iquadfams=1:length(quadfams)
  KQ=cellfun(@(a)a.PolynomB(2),arc(quadfams{iquadfams}),'un',1);
  variabs=[variabs, atVariableBuilder(arc,...
    {@(r,DKquad)setcellstruct(r,'PolynomB',quadfams{iquadfams},KQ+DKquad,1,2)},...
    {[1e-8]})]; %#ok<*AGROW>
end

ConstrQX=struct(...
    'Fun',@(~,ld,~)mux(ld),...
    'Weight',1,...
    'RefPoints',[1:length(arc)+1],...
    'Min',tune(1),...
    'Max',tune(1));

ConstrQY=struct(...
    'Fun',@(~,ld,~)muy(ld),...
    'Weight',1,...
    'RefPoints',[1:length(arc)+1],...
    'Min',tune(2),...
    'Max',tune(2));

% tol=1e-6;
% arctune0=atmatch(arc,variabs,[ConstrQX ConstrQY],tol,5,3);%,@lsqnonlin); %

 tol=1e-8;
 arctune0=arc;
 arctune0=atmatch(arctune0,variabs,[ConstrQX ConstrQY],tol,10,3,@lsqnonlin); %);%
 arctune0=atmatch(arctune0,variabs,[ConstrQX ConstrQY],tol,50,3); %);%

return

function m=muy(lindata)

m=lindata(end).mu(2)/2/pi;

return

function m=mux(lindata)

m=lindata(end).mu(1)/2/pi;

return
