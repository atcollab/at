function plotdata = plotRDT(lindata,ring,dpp, varargin)
%
%  plotRDT plots the absolute value of the hamiltonian terms
%  plotRDT must be used with atplot:
%  
%  atplot(ring,@plotRDT,'geometric1') plots the first order geometric terms
%  atplot(ring,@plotRDT,'chromatic') plots the chromatic terms
%  atplot(ring,@plotRDT,'geometric2') plots the second order geometric terms
%  atplot(ring,@plotRDT,'coupling') plots the coupling terms
%  
%  see also: computeRDT, atplot

naddvar=length(varargin);
chromatic=0;
coupling=0;
geometric1=0;
geometric2=0;
tuneshifts=0;
OnlyFirst=0;
if(naddvar>0)
    %for ii=1:naddvar
    switch varargin{1}
        case 'chromatic'
            chromatic=1;
        case 'coupling'
            coupling=1;
        case 'geometric1'
            geometric1=1;
        case 'geometric2'
            geometric2=1;
        case 'tuneshifts'
            tuneshifts=1;
        otherwise
            disp('The first input must be one of these:');
            disp('''chromatic'', ''coupling'', ''geometric1'', ''geometric2''');
            disp('your input will be considered ''geometric1''');
            geometric1=1;
    end
    if naddvar>1
        switch varargin{2}
            case 'OnlyFirstOrders'
                OnlyFirst=1;
            otherwise
                disp(['The second input can be only:']);
                disp('''OnlyFirstOrders''');
                disp('your second input will be ignored');
        end
    end
else
    chromatic=0;
    coupling=0;
    geometric1=1;
    geometric2=0;
    tuneshifts=0;
    OnlyFirst=0;
end

indSext=findcells(ring,'Class','Sextupole');
indSext=unique(sort([ 1, length(ring),indSext+1]));
indOct=findcells(ring,'Class','Octupole');
indOct=unique(sort([1,length(ring),indOct+1]));
indSO=unique(sort([indSext,indOct]));
indQuad=findcells(ring,'Class','Quadrupole');
indQuad=unique(sort([1,length(ring),indQuad+1]));
indDip=findcells(ring,'Class','Bend');
indDip=unique(sort([1,length(ring),indDip+1]));
indTot=unique(sort([ indSext, indOct, indQuad, indDip ]));
if(geometric1 || geometric2 || tuneshifts)
    if(geometric1)
        RDT=computeRDT(ring,indSO,'geometric1');
    else
        if OnlyFirst==1
            RDT=computeRDT(ring,indOct,varargin{1},'OnlyFirstOrders');
        else
            RDT=computeRDT(ring,indSO,varargin{1});
        end
    end
else
    RDT=computeRDT(ring,indTot,varargin{1});
end

Spos=findspos(ring,1:length(ring));
if(chromatic)
    h11001=[]; h00111=[]; h20001=[]; h00201=[]; h10002=[];
    h11001p=abs(arrayfun(@(x)x.h11001,RDT));
    h00111p=abs(arrayfun(@(x)x.h00111,RDT));
    h20001p=abs(arrayfun(@(x)x.h20001,RDT));
    h00201p=abs(arrayfun(@(x)x.h00201,RDT));
    h10002p=abs(arrayfun(@(x)x.h10002,RDT));
    for(ii=1:length(indTot)-1)
        h11001=[h11001; repmat(h11001p(ii),indTot(ii+1)-indTot(ii),1)];
        h00111=[h00111; repmat(h00111p(ii),indTot(ii+1)-indTot(ii),1)];
        h20001=[h20001; repmat(h20001p(ii),indTot(ii+1)-indTot(ii),1)];
        h00201=[h00201; repmat(h00201p(ii),indTot(ii+1)-indTot(ii),1)];
        h10002=[h10002; repmat(h10002p(ii),indTot(ii+1)-indTot(ii),1)];
    end
    h11001=[h11001; repmat(h11001p(length(indTot)),2,1)];
    h00111=[h00111; repmat(h00111p(length(indTot)),2,1)];
    h20001=[h20001; repmat(h20001p(length(indTot)),2,1)];
    h00201=[h00201; repmat(h00201p(length(indTot)),2,1)];
    h10002=[h10002; repmat(h10002p(length(indTot)),2,1)];
    
    plotdata(1).values=[ h11001 h00111 h20001 h00201 ];
    plotdata(2).values=[ h10002 ];
    plotdata(1).labels={ 'h11001','h00111','h20001','h00201' };
    plotdata(2).labels={'h10002'};
    plotdata(1).axislabel='absolute value of RDTs';
    plotdata(2).axislabel='absolute value of h10002';
end
if(coupling)
    h10010=[]; h10100=[];
    h10010p=abs(arrayfun(@(x)x.h10010,RDT));
    h10100p=abs(arrayfun(@(x)x.h10100,RDT));
    for(ii=1:length(indTot)-1)
        h10010=[h10010; repmat(h10010p(ii),indTot(ii+1)-indTot(ii),1)];
        h10100=[h10100; repmat(h10100p(ii),indTot(ii+1)-indTot(ii),1)];
    end
    h10010=[h10010; repmat(h10010p(length(indTot)),2,1)];
    h10100=[h10100; repmat(h10100p(length(indTot)),2,1)];
    
    plotdata(1).values=[h10010 h10100];
    plotdata(1).labels={'h10010','h10100'};
    plotdata(1).axislabel='absolute value of RDTs';
end
if(geometric1)
    h21000=[]; h30000=[]; h10110=[]; h10020=[]; h10200=[];
    h21000p=abs(arrayfun(@(x)x.h21000,RDT));
    h30000p=abs(arrayfun(@(x)x.h30000,RDT));
    h10110p=abs(arrayfun(@(x)x.h10110,RDT));
    h10020p=abs(arrayfun(@(x)x.h10020,RDT));
    h10200p=abs(arrayfun(@(x)x.h10200,RDT));
    for(ii=1:length(indSO)-1)
        h21000=[h21000; repmat(h21000p(ii),indSO(ii+1)-indSO(ii),1)];
        h30000=[h30000; repmat(h30000p(ii),indSO(ii+1)-indSO(ii),1)];
        h10110=[h10110; repmat(h10110p(ii),indSO(ii+1)-indSO(ii),1)];
        h10020=[h10020; repmat(h10020p(ii),indSO(ii+1)-indSO(ii),1)];
        h10200=[h10200; repmat(h10200p(ii),indSO(ii+1)-indSO(ii),1)];
    end
    h21000=[h21000; repmat(h21000p(ii),2,1)];
    h30000=[h30000; repmat(h30000p(ii),2,1)];
    h10110=[h10110; repmat(h10110p(ii),2,1)];
    h10020=[h10020; repmat(h10020p(ii),2,1)];
    h10200=[h10200; repmat(h10200p(ii),2,1)];
    plotdata(1).values=[h21000 h30000 h10110 h10020 h10200];
    plotdata(1).labels={'h21000','h30000','h10110','h10020','h10200'};
    plotdata(1).axislabel='absolute value of RDTs';
end
if(geometric2)
    h22000=[]; h11110=[]; h00220=[]; h31000=[]; h40000=[]; h20110=[]; h11200=[]; h20020=[]; h20200=[]; h00310=[]; h00400=[];
    h22000p=abs(arrayfun(@(x)x.h22000,RDT));
    h11110p=abs(arrayfun(@(x)x.h11110,RDT));
    h00220p=abs(arrayfun(@(x)x.h00220,RDT));
    h31000p=abs(arrayfun(@(x)x.h31000,RDT));
    h40000p=abs(arrayfun(@(x)x.h40000,RDT));
    h20110p=abs(arrayfun(@(x)x.h20110,RDT));
    h11200p=abs(arrayfun(@(x)x.h11200,RDT));
    h20020p=abs(arrayfun(@(x)x.h20020,RDT));
    h20200p=abs(arrayfun(@(x)x.h20200,RDT));
    h00310p=abs(arrayfun(@(x)x.h00310,RDT));
    h00400p=abs(arrayfun(@(x)x.h00400,RDT));
    if (~OnlyFirst)
    for(ii=1:length(indSO)-1)
        h22000=[h22000; repmat(h22000p(ii),indSO(ii+1)-indSO(ii),1)];
        h11110=[h11110; repmat(h11110p(ii),indSO(ii+1)-indSO(ii),1)];
        h00220=[h00220; repmat(h00220p(ii),indSO(ii+1)-indSO(ii),1)];
        h31000=[h31000; repmat(h31000p(ii),indSO(ii+1)-indSO(ii),1)];
        h40000=[h40000; repmat(h40000p(ii),indSO(ii+1)-indSO(ii),1)];
        h20110=[h20110; repmat(h20110p(ii),indSO(ii+1)-indSO(ii),1)];
        h11200=[h11200; repmat(h11200p(ii),indSO(ii+1)-indSO(ii),1)];
        h20020=[h20020; repmat(h20020p(ii),indSO(ii+1)-indSO(ii),1)];
        h20200=[h20200; repmat(h20200p(ii),indSO(ii+1)-indSO(ii),1)];
        h00310=[h00310; repmat(h00310p(ii),indSO(ii+1)-indSO(ii),1)];
        h00400=[h00400; repmat(h00400p(ii),indSO(ii+1)-indSO(ii),1)];
    end
    h22000=[h22000; repmat(h22000p(ii),2,1)];
    h11110=[h11110; repmat(h11110p(ii),2,1)];
    h00220=[h00220; repmat(h00220p(ii),2,1)];
    h31000=[h31000; repmat(h31000p(ii),2,1)];
    h40000=[h40000; repmat(h40000p(ii),2,1)];
    h20110=[h20110; repmat(h20110p(ii),2,1)];
    h11200=[h11200; repmat(h11200p(ii),2,1)];
    h20020=[h20020; repmat(h20020p(ii),2,1)];
    h20200=[h20200; repmat(h20200p(ii),2,1)];
    h00310=[h00310; repmat(h00310p(ii),2,1)];
    h00400=[h00400; repmat(h00400p(ii),2,1)];
    else
        for(ii=1:length(indOct)-1)
        h22000=[h22000; repmat(h22000p(ii),indOct(ii+1)-indOct(ii),1)];
        h11110=[h11110; repmat(h11110p(ii),indOct(ii+1)-indOct(ii),1)];
        h00220=[h00220; repmat(h00220p(ii),indOct(ii+1)-indOct(ii),1)];
        h31000=[h31000; repmat(h31000p(ii),indOct(ii+1)-indOct(ii),1)];
        h40000=[h40000; repmat(h40000p(ii),indOct(ii+1)-indOct(ii),1)];
        h20110=[h20110; repmat(h20110p(ii),indOct(ii+1)-indOct(ii),1)];
        h11200=[h11200; repmat(h11200p(ii),indOct(ii+1)-indOct(ii),1)];
        h20020=[h20020; repmat(h20020p(ii),indOct(ii+1)-indOct(ii),1)];
        h20200=[h20200; repmat(h20200p(ii),indOct(ii+1)-indOct(ii),1)];
        h00310=[h00310; repmat(h00310p(ii),indOct(ii+1)-indOct(ii),1)];
        h00400=[h00400; repmat(h00400p(ii),indOct(ii+1)-indOct(ii),1)];
    end
    h22000=[h22000; repmat(h22000p(ii),2,1)];
    h11110=[h11110; repmat(h11110p(ii),2,1)];
    h00220=[h00220; repmat(h00220p(ii),2,1)];
    h31000=[h31000; repmat(h31000p(ii),2,1)];
    h40000=[h40000; repmat(h40000p(ii),2,1)];
    h20110=[h20110; repmat(h20110p(ii),2,1)];
    h11200=[h11200; repmat(h11200p(ii),2,1)];
    h20020=[h20020; repmat(h20020p(ii),2,1)];
    h20200=[h20200; repmat(h20200p(ii),2,1)];
    h00310=[h00310; repmat(h00310p(ii),2,1)];
    h00400=[h00400; repmat(h00400p(ii),2,1)];
    end
    plotdata(1).values=[h22000 h11110 h00220 h31000 h40000 h20110 h11200 h20020 h20200 h00310 h00400];
    plotdata(1).labels={'h22000','h11110','h00220','h31000','h40000','h20110','h11200','h20020','h20200','h00310','h00400'};
    plotdata(1).axislabel='absolute value of RDTs';
end



