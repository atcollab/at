function lmt=CreateLatexMagnetsTableCompSingleClass(r0,r1,atclass,tabname)
% generate a latex input file for a table of magnet fields, names and length
% lmt is a string containing the same information.
%
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

switch atclass
    case {'Dipole','Bend'}
        ind=1;
    case {'Quadrupole'}
        ind=2;
    case {'Sextupole'}
        ind=3;
end    

tabhead=['\\begin{table}[!h]\n'...
    '\\caption{' strrep(tabname,'_',' ') '}\n'...
    '\\begin{center}\n'...
    '\\begin{tabular}{rcc}\n'...
    '\\toprule\n'...
    ' & before &after \\\\\n'...
    '\\midrule\n'...
    'name &$b_' num2str(ind) '\\,[T/m^' num2str(ind-1) ']$ &$b_' num2str(ind) '\\,[T/m^' num2str(ind-1) ']$   \\\\\n'...
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
%E0=r{1}.Energy;

E0=6.04e9;

Brho=3.3356*E0/1e9;

b0=findcells(r0,'Class',atclass);
b1=findcells(r1,'Class',atclass);


% get names, b1,b2,betax,etax
% names=cellfun(@(x)x.FamName,r(mag),'un',0);
% [~,i]=unique(names);
% [mag,im]=sort(mag(i));
% mag=mag(end:-1:1);
% colors=colors(im(end:-1:1),:);
% 

L=cellfun(@(x)x.Length,r0(b0));

names=cellfun(@(x)x.FamName,r0(b0),'un',0);

B3_0=cellfun(@(x)x.PolynomB(ind),r0(b0),'un',0);
B3_1=cellfun(@(x)x.PolynomB(ind),r1(b1),'un',0);

B3_0=cell2mat(B3_0)*Brho;
B3_1=cell2mat(B3_1)*Brho;

% print to string
leftparstr=repmat('{',size(L));
rightparstr=repmat('}',size(L));
andstr=repmat('&',size(L));
nlstr=repmat('\\\\ \n',size(L));

namesl=zeros(length(names),10);
for i=1:length(names)
    LN=length(names{i});
    if LN<10
        namesl(i,:)=[names{i} repmat(' ',1,(10-LN))];
    end
end

lmt=[namesl andstr ...
     num2str(B3_0,'%3.1f') andstr...
     num2str(B3_1,'%3.1f') nlstr];

lmt=reshape(lmt',1,numel(lmt));

lmt=strrep(lmt,'0.000}','    }');
lmt=strrep(lmt,'0.00}','    }');
lmt=strrep(lmt,'0.0}','    }');

% print to file
fff=fopen([tabname '.tex'],'w+');

fprintf(fff,tabhead);
fprintf(fff,lmt');
fprintf(fff,tabtail);

fclose('all')

close all
clear all

return