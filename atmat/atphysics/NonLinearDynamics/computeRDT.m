function RDT=computeRDT(ring, index, varargin)
%COMPUTERDT Computes Hamiltonian resonance driving terms (RDTs)
%   This function calls RDTElegantAT mex function and returns the
%   hamiltonian resonance driving terms, using the elegant c++ 
%   function computeDrivingTerms().
%   
%   RDT=computeRDT(ring, index, varargin)
%
%   ring is the AT lattice
%   index is the vector of indexes where one wants to compute RDTs
%   The additional arguments can be up to five strings:
%   chromatic, coupling, geometric1, geometric2 and tuneshifts
%   
%   example:
%   RDT=computeRDT(ring, indexBPM, 'geometric1', 'tuneshifts');
%   creates an array of structs (the length of the array is the number of 
%   indexes where you want to compute driving terms) with first order
%   geometric driving terms and tune shifts with amplitude.
%   The driving terms are complex numbers, the tune shifts are real.
%

naddvar=length(varargin);
chromatic=0;
coupling=0;
geometric1=0;
geometric2=0;
tuneshifts=0;
for ii=1:naddvar
    switch varargin{ii}
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
            disp(['The input number ' num2str(ii+2) ' must be one of these:']);
            disp('''chromatic'', ''coupling'', ''geometric1'',''geometric2'', ''tuneshifts''');
            disp('your input will be ignored');
    end
end

if naddvar==0
    chromatic=1;
    coupling=1;
    geometric1=1;
    geometric2=1;
    tuneshifts=1;
end


indDQSO=findcells(ring,'Class','Bend','Quadrupole','Sextupole','Octupole','Multipole');

[~,AVEBETA,AVEMU,AVEDISP,~,~]=atavedata(ring,0,1:length(ring));

Lin=atlinopt(ring,0,1:length(ring));

%create input arguments for the mex function
sIndex=findspos(ring,indDQSO);
s=findspos(ring,1:length(ring));
sEnd=findspos(ring,length(ring)+1);
betax=AVEBETA(indDQSO,1);
betay=AVEBETA(indDQSO,2);
etax=AVEDISP(indDQSO,1);
phix=AVEMU(indDQSO,1);
phiy=AVEMU(indDQSO,2);
a2L=getcellstruct(ring,'PolynomA',indDQSO,1,2).*getcellstruct(ring,'Length',indDQSO);
a2L(isnan(a2L))=0;
b2L=getcellstruct(ring,'PolynomB',indDQSO,1,2).*getcellstruct(ring,'Length',indDQSO);
b2L(isnan(b2L))=0;
b3L=getcellstruct(ring,'PolynomB',indDQSO,1,3).*getcellstruct(ring,'Length',indDQSO);
b3L(isnan(b3L))=0;
b4L=getcellstruct(ring,'PolynomB',indDQSO,1,4).*getcellstruct(ring,'Length',indDQSO);
b4L(isnan(b4L))=0;
Mux=Lin(length(ring)).mu(1);
Tunex=Mux/2/pi;
Muy=Lin(length(ring)).mu(2);
Tuney=Muy/2/pi;
nElem=length(indDQSO);

for ii=1:length(index)
    FromindexDQSO=sum(indDQSO<index(ii))+1;
    betax_Fromindex=[betax(FromindexDQSO:end);betax(1:FromindexDQSO-1)];
    betay_Fromindex=[betay(FromindexDQSO:end);betay(1:FromindexDQSO-1)];
    etax_Fromindex=[etax(FromindexDQSO:end);etax(1:FromindexDQSO-1)];
    phix_Fromindex=[phix(FromindexDQSO:end)-AVEMU(index(ii),1);phix(1:FromindexDQSO-1)+Mux-AVEMU(index(ii),1)];
    phiy_Fromindex=[phiy(FromindexDQSO:end)-AVEMU(index(ii),2);phiy(1:FromindexDQSO-1)+Muy-AVEMU(index(ii),2)];
    s_Fromindex=[sIndex(FromindexDQSO:end)-s(index(ii)),sIndex(1:FromindexDQSO-1)+sEnd-s(index(ii))];
    a2L_Fromindex=[a2L(FromindexDQSO:end);a2L(1:FromindexDQSO-1)];
    b2L_Fromindex=[b2L(FromindexDQSO:end);b2L(1:FromindexDQSO-1)];
    b3L_Fromindex=[b3L(FromindexDQSO:end);b3L(1:FromindexDQSO-1)];
    b4L_Fromindex=[b4L(FromindexDQSO:end);b4L(1:FromindexDQSO-1)];
        [ReRDT, ImRDT, TSwA]=RDTelegantAT(s_Fromindex,betax_Fromindex,betay_Fromindex,...
        etax_Fromindex,phix_Fromindex,phiy_Fromindex,a2L_Fromindex,b2L_Fromindex,...
        b3L_Fromindex,b4L_Fromindex,Tunex,Tuney,nElem,...
        chromatic,coupling,geometric1,geometric2,tuneshifts);
    %chromatic
    if(chromatic)
    RDT(ii).h11001=ReRDT(6)+1i*ImRDT(6);
    RDT(ii).h00111=ReRDT(7)+1i*ImRDT(7);
    RDT(ii).h20001=ReRDT(8)+1i*ImRDT(8);
    RDT(ii).h00201=ReRDT(9)+1i*ImRDT(9);
    RDT(ii).h10002=ReRDT(10)+1i*ImRDT(10);
    end
    %coupling
    if(coupling)
    RDT(ii).h10010=ReRDT(11)+1i*ImRDT(11);
	RDT(ii).h10100=ReRDT(12)+1i*ImRDT(12);
    end
    %geometric1
    if(geometric1)
	RDT(ii).h21000=ReRDT(1)+1i*ImRDT(1);
	RDT(ii).h30000=ReRDT(2)+1i*ImRDT(2);
	RDT(ii).h10110=ReRDT(3)+1i*ImRDT(3);
	RDT(ii).h10020=ReRDT(4)+1i*ImRDT(4);
	RDT(ii).h10200=ReRDT(5)+1i*ImRDT(5);
    end
    %geometric2
    if(geometric2)
    RDT(ii).h22000=ReRDT(13)+1i*ImRDT(13);
	RDT(ii).h11110=ReRDT(14)+1i*ImRDT(14);
	RDT(ii).h00220=ReRDT(15)+1i*ImRDT(15);
	RDT(ii).h31000=ReRDT(16)+1i*ImRDT(16);
	RDT(ii).h40000=ReRDT(17)+1i*ImRDT(17);
	RDT(ii).h20110=ReRDT(18)+1i*ImRDT(18);
	RDT(ii).h11200=ReRDT(19)+1i*ImRDT(19);
	RDT(ii).h20020=ReRDT(20)+1i*ImRDT(20);
	RDT(ii).h20200=ReRDT(21)+1i*ImRDT(21);
	RDT(ii).h00310=ReRDT(22)+1i*ImRDT(22);
	RDT(ii).h00400=ReRDT(23)+1i*ImRDT(23);
    end
    %tuneshifts
    if(tuneshifts)
    RDT(ii).dnux_dJx=TSwA(1);
	RDT(ii).dnux_dJy=TSwA(2);
	RDT(ii).dnuy_dJy=TSwA(3);
    end
end