function lmt=CreateLatexMagnetsTable(r,tabname,E0)
% function lmt=CreateLatexMagnetsTable(r,tabname)
%
% r= at lattice
% tabname= name of table out file. (label will be refto_tabname)
% 
% generate a latex input file for a table of magnet fields, names and length
% lmt is a string containing the same information.
%
% needs package
% \usepackage{xcolor}
% 
% to hide columns: 
% \newcolumntype{H}{>{\setbox0=\hbox\bgroup}c<{\egroup}@{}}
% 
% in latex the table may be inserted as
% 
% \input{tabname}
% 

%bety
%colori come in synoptic per i magneti.
%theta, R0, B0
%vuoto invece di zero
%solo first half cell

% r=mergedoubles(r,'QF0');
% %r=mergedoubles(r,'BPI2');
% indx=findcells(r,'FamName','BPI2');
% r0=r;
% r=setcellstruct(r,'Length',indx,getcellstruct(r0,'Length',indx)*2); %empty
% r=setcellstruct(r,'BendingAngle',indx,getcellstruct(r0,'BendingAngle',indx)*2); %empty
% r(indx([2:2:end]))=[];
% 
% indx=findcells(r,'FamName','BPI3');
% r0=r;
% r=setcellstruct(r,'Length',indx,getcellstruct(r0,'Length',indx)*2); %empty
% r=setcellstruct(r,'BendingAngle',indx,getcellstruct(r0,'BendingAngle',indx)*2); %empty
% r(indx([2:2:end]))=[];
% 
%r=mergedoubles(r,'BPI3');
%r=mergedoubles(r,'SF1A');
%r=mergedoubles(r,'SD1A');

[l,~,~]=atlinopt(r,0,1:length(r));


tabhead=['\\begin{table}[!h]\n'...
    '\\caption{' tabname '}\n'...
    '\\begin{center}\n'...
    '\\begin{tabular}{rcccccccccc}\n'...
    '\\toprule\n'...
    'name & L [m] & $\\rho\\,[m]$ & $\\theta\\,[rad]$& $B_0\\,[T]$ &$b_2\\,[T/m]$ &$b_3\\,[T/m^2]$ &$b_4\\,[T/m^3]$ & $\\beta_x\\,[m]$& $\\beta_y\\,[m]$ &$\\eta_x\\,[m]$ \\\\\n'...
    '\\midrule\n'...
    ];

tabtail=[...
    '\\bottomrule\n'...
    '\\end{tabular}\n'....
    '\\end{center}\n'...
    '\\label{refto_' tabname '}\n'...
    '\\end{table}\n'...
    ...
    ];

% get Brho
if nargin <3
    try
        E0=r{1}.Energy;
    catch
        try
            E0=GLOBVAL.E0;
        catch
            warning('no energy defined set  E0=0.')
            E0=0;
        end
    end
end

Brho=3.3356*E0/1e9;

if isfield(r{1},'Class')
    b=findcells(r,'Class','Bend');
    q=findcells(r,'Class','Quadrupole');
    s=findcells(r,'Class','Sextupole');
    o=findcells(r,'Class','Octupole');
else
    b=findcells(r,'BetaCode','DI');
    q=findcells(r,'BetaCode','QP');
    s=findcells(r,'BetaCode','SX');
    o=findcells(r,'BetaCode','OC');
end

% only half cell
b=b(b<ceil(length(r)/2+1));
q=q(q<ceil(length(r)/2+1));
s=s(s<ceil(length(r)/2+1));
o=o(o<ceil(length(r)/2+1));

% get magnets
colors=[ repmat('\\color{blue!80!white}  ',length(b),1);...
         repmat('\\color{red!80!black}   ',length(q),1);...
         repmat('\\color{green!80!black} ',length(s),1);...
         repmat('\\color{yellow!60!black}',length(o),1)];

[mag,magord]=sort([b,q,s,o]);

colors=colors(magord,:);


% get names, b1,b2,betax,etax
% names=cellfun(@(x)x.FamName,r(mag),'un',0);
% [~,i]=unique(names);
% [mag,im]=sort(mag(i));
% mag=mag(end:-1:1);
% colors=colors(im(end:-1:1),:);
% 

dips=findcells(r(mag),'BendingAngle');

L=cellfun(@(x)x.Length,r(mag));

names=cellfun(@(x)x.FamName,r(mag),'un',0);
B0=cellfun(@(x)x.PolynomB(1),r(mag),'un',0);
B0(dips)=cellfun(@(x)x.BendingAngle/x.Length,r(mag(dips)),'un',0);
B1=cellfun(@(x)x.PolynomB(2),r(mag),'un',0);
B2=cellfun(@(x)x.PolynomB(3),r(mag),'un',0);
B3_=cellfun(@(x)x.PolynomB(4),r(o),'un',0);

betaX=arrayfun(@(x)x.beta(1),l(mag),'un',0)';
betaY=arrayfun(@(x)x.beta(2),l(mag),'un',0)';
disp=arrayfun(@(x)x.Dispersion(1),l(mag),'un',0)';

%names=cell2mat(names);
%L=cell2mat(L);

TH=cell2mat(B0);
B0=cell2mat(B0)*Brho;
B1=cell2mat(B1)*Brho;
B2=cell2mat(B2)*Brho;
B3_=cell2mat(B3_)*Brho;
B3=zeros(size(B2));
B3(mag==o)=B3_;

betaX=cell2mat(betaX);
betaY=cell2mat(betaY);
disp=cell2mat(disp);

RH=L./TH;
RH(isinf(RH))=0;

% print to string
leftparstr=repmat('{',size(B0));
rightparstr=repmat('}',size(B0));
andstr=repmat('&',size(B0));
nlstr=repmat('\\\\ \n',size(B0));


namesl=zeros(length(names),10);
for i=1:length(names)
    LN=length(names{i});
    if LN<10
        namesl(i,:)=[names{i} repmat(' ',1,(10-LN))];
    end
end

lmt=[...
    leftparstr colors namesl rightparstr andstr ...
    num2str(L,'%1.3f') andstr...
    leftparstr colors num2str(RH,'%1.3f') rightparstr andstr...
    leftparstr colors num2str(TH,'%1.3f') rightparstr andstr...
    leftparstr colors num2str(B0,'%1.3f') rightparstr andstr...
    leftparstr colors num2str(B1,'%2.1f') rightparstr andstr...
    leftparstr colors num2str(B2,'%3.1f') rightparstr andstr...
    leftparstr colors num2str(B3,'%3.1f') rightparstr andstr...
    num2str(betaX,'%2.2f') andstr...
    num2str(betaY,'%2.2f') andstr...
    num2str(disp,'%1.3f') nlstr];

lmt=reshape(lmt',1,numel(lmt));

lmt=strrep(lmt,'0.000}','    }');
%lmt=strrep(lmt,'0.00}','    }');
lmt=strrep(lmt,'}  0.0}','}  }');
lmt=strrep(lmt,'}   0.0}','}   }');
lmt=strrep(lmt,'}    0.0}','}   }');
lmt=strrep(lmt,'}     0.0}','}   }');
lmt=strrep(lmt,'}      0.0}','}   }');
lmt=strrep(lmt,'}       0.0}','}   }');
lmt=strrep(lmt,'}        0.0}','}   }');

% print to file
fff=fopen([tabname '.tex'],'w+');

fprintf(fff,tabhead);
fprintf(fff,lmt');
fprintf(fff,tabtail);

fclose('all')

close all
clear all

return