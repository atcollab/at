function ModelRM...
    =getresponsematrices(...
    r0,...          %1 AT lattice
    indBPM,...      %2 bpm indexes in at lattice
    indHorCor,...   %3 h cor indexes
    indVerCor,...   %4 v cor indexes
    indSkewCor,...  %5 skew cor indexes
    indQuadCor,...  %6 quad cor indexes
    indSextCor,...  %7 sext cor indexes
    inCOD,...       %7 initial coordinates
    rmsel...        %8 specifiy rm to be computed
    )
% function ModelRM...
%     =getorbitbetaphasedispersionresponse(...
%     r0,...          %1 AT lattice
%     indBPM,...      %2 bpm indexes in at lattice
%     indHorCor,...   %3 h cor indexes
%     indVerCor,...   %4 v cor indexes
%     indSkewCor,...  %5 skew cor indexes
%     indQuadCor,...  %6 quad cor indexes
%     indSextCor,...  %7 sext cor indexes
%     inCOD,...)      %7 [0 0 0 0 0 0]' inCOD
%     rmsel)          %8 [1 2 ... 12] rm to compute (default: all)
%
% computes lattice Response Matrix for :
%
% rmsel                                   output structure is:
%
%  1     D  orbit (6D) / D h steerer       ModelRM.OrbHCor
%  2     D  orbit (6D) / D v steerer       ModelRM.OrbHCor
%  3     D  orbit (6D) / D dpp             ModelRM.OrbHDpp
%        D  orbit (6D) / D dpp             ModelRM.OrbVDpp
%  4     D  trajectory (6D) / D h steerer  ModelRM.TrajHCor
%  5     D  trajectory (6D) / D v steerer  ModelRM.TrajVCor
%  6     D  trajectory (6D) / D dpp        ModelRM.TrajHDpp
%        D  trajectory (6D) / D dpp        ModelRM.TrajVDpp
%  7     D  dispersion (6D) / D h steerer  ModelRM.DispHCor
%  8     D  dispersion (6D) / D v steerer  ModelRM.DispVCor
%  9     D  dispersion (6D) / D dpp        ModelRM.DispHDpp
%        D  dispersion (6D) / D dpp        ModelRM.DispVDpp
% 10     D  dispersion (6D) / D norm quad  ModelRM.DispQCor
% 11     D  dispersion (6D) / D skew quad  ModelRM.DispSCor
% 12     D  tune            / D norm quad  ModelRM.TuneQCor
%
% for the dispersion measurement the RF freq. is changed (finddispersion6Err).
%
%see also: findrespmat findorbit6Err findtrajectory6Err finddispersion6Err

%[GRM,GDRM]=getGirderRM(r0,indBPM);
%disp(' --- computed Girders RM --- ')

kval=1e-4;
delta=1e-3;

if nargin<9
    rmsel=1:12;
end

if nargin<8
    inCOD=[0 0 0 0 0 0]';
end

% orbit RM dpp
alpha=mcf(r0);
indrfc=find(atgetcells(r0,'Frequency'));
f0=r0{indrfc(1)}.Frequency;

for ir=1:length(rmsel)
    switch rmsel(ir)
        case 1
            % orbit RM
            ormH=findrespmat(r0,indBPM,indHorCor,kval,'PolynomB',1,1,'findorbit6Err',inCOD);
            ormH{1}=ormH{1}./kval;
            ormH{2}=ormH{2}./kval;
            ormH{3}=ormH{3}./kval;
            ormH{4}=ormH{4}./kval;
            disp(' --- computed H orbit RM from model --- ');
            % store data
            ModelRM.OrbHCor=ormH;
        case 2
            
            ormV=findrespmat(r0,indBPM,indVerCor,kval,'PolynomA',1,1,'findorbit6Err',inCOD);
            ormV{1}=ormV{1}./kval;
            ormV{2}=ormV{2}./kval;
            ormV{3}=ormV{3}./kval;
            ormV{4}=ormV{4}./kval;
            disp(' --- computed V orbit RM from model --- ');
            
            % store data
            ModelRM.OrbVCor=ormV;
            
        case 3
            % plus delta
            RINGp=atsetfieldvalues(r0,indrfc,'Frequency',f0-alpha*(+delta)*f0);
            o=findorbit6Err(RINGp,indBPM,inCOD);
            oxPdpp=o(1,:);
            oyPdpp=o(3,:);
            
            RINGm=atsetfieldvalues(r0,indrfc,'Frequency',f0-alpha*(-delta)*f0);
            o=findorbit6Err(RINGm,indBPM,inCOD);
            oxMdpp=o(1,:);
            oyMdpp=o(3,:);
            
            dppH=(oxPdpp-oxMdpp)/2/delta;
            dppV=(oyPdpp-oyMdpp)/2/delta;
            
            disp(' --- computed orbit/dpp RM from model --- ');
            
            % store data
            
            ModelRM.OrbHDPP=dppH;
            ModelRM.OrbVDPP=dppV;
            
        case 4
            % trajectory
            RMTH=findrespmat(r0,indBPM,indHorCor,kval,'PolynomB',1,1,'findtrajectory6Err',inCOD);
            RMTH{1}=RMTH{1}./kval;
            RMTH{2}=RMTH{2}./kval;
            RMTH{3}=RMTH{3}./kval;
            RMTH{4}=RMTH{4}./kval;
            disp(' --- computed horizontal trajectory RM from model --- ');
            % store data
            ModelRM.TrajHCor=RMTH;
            
        case 5
            
            RMTV=findrespmat(r0,indBPM,indVerCor,kval,'PolynomA',1,1,'findtrajectory6Err',inCOD);
            RMTV{1}=RMTV{1}./kval;
            RMTV{2}=RMTV{2}./kval;
            RMTV{3}=RMTV{3}./kval;
            RMTV{4}=RMTV{4}./kval;
            disp(' --- computed vertical trajectory RM from model --- ');
            
            % store data
            ModelRM.TrajVCor=RMTV;
            
            
        case 6
            % orbit RM dpp
            inCODPdpp=inCOD;
            inCODPdpp(5)=delta;% linepass is used.
            o=findtrajectory6Err(r0,indBPM,inCODPdpp);
            oxPdpp=o(1,:);
            oyPdpp=o(3,:);
            inCODMdpp=inCOD;
            inCODMdpp(5)=-delta;
            o=findtrajectory6Err(r0,indBPM,inCODMdpp);
            oxMdpp=o(1,:);
            oyMdpp=o(3,:);
            dppH=(oxPdpp-oxMdpp)/2/delta;
            dppV=(oyPdpp-oyMdpp)/2/delta;
            disp(' --- computed trajectory/dpp RM from model --- ');
            
            % store data
            ModelRM.TrajHDPP=dppH;
            ModelRM.TrajVDPP=dppV;
            
        case 7
            
            % dispersion RM steerers
            
            drmH=findrespmat(r0,indBPM,indHorCor,kval,'PolynomB',1,1,'finddispersion6Err',indrfc,alpha,delta,inCOD);
            drmH{1}=drmH{1}./kval;
            drmH{2}=drmH{2}./kval;
            drmH{3}=drmH{3}./kval;
            drmH{4}=drmH{4}./kval;
            disp(' --- computed H dispersion RM steerer h from model --- ');
            % store data
            ModelRM.DispHCor=drmH;
            
        case 8
            
            
            drmV=findrespmat(r0,indBPM,indHorCor,kval,'PolynomA',1,1,'finddispersion6Err',indrfc,alpha,delta,inCOD);
            drmV{1}=drmV{1}./kval;
            drmV{2}=drmV{2}./kval;
            drmV{3}=drmV{3}./kval;
            drmV{4}=drmV{4}./kval;
            disp(' --- computed H dispersion RM steerer v from model --- ');
            
            % store data
            ModelRM.DispVCor=drmV;
        case 9
            
            % plus delta
            RINGp=atsetfieldvalues(r0,indrfc,'Frequency',f0-alpha*(+delta)*f0);
            o=finddispersion6Err(RINGp,indBPM,indrfc,alpha,delta,inCOD);
            oxPdpp=o(1,:);
            oyPdpp=o(3,:);
            
            RINGm=atsetfieldvalues(r0,indrfc,'Frequency',f0-alpha*(-delta)*f0);
            o=finddispersion6Err(RINGm,indBPM,indrfc,alpha,delta,inCOD);
            oxMdpp=o(1,:);
            oyMdpp=o(3,:);
            
            dppH=(oxPdpp-oxMdpp)/2/delta;
            dppV=(oyPdpp-oyMdpp)/2/delta;
            disp(' --- computed dispersion/dpp RM from model --- ');
            
            % store data
            ModelRM.DispHDPP=dppH;
            ModelRM.DispVDPP=dppV;
            
        case 10
            % dispersion RM quadrupoles
            
            drmQ=findrespmat(r0,indBPM,indQuadCor,kval,'PolynomB',1,2,'finddispersion6Err',indrfc,alpha,delta,inCOD);
            drmQ{1}=drmQ{1}./kval;
            drmQ{2}=drmQ{2}./kval;
            drmQ{3}=drmQ{3}./kval;
            drmQ{4}=drmQ{4}./kval;
            disp(' --- computed H dispersion RM norm quad from model --- ');
            % store data
            ModelRM.DispQCor=drmQ;
            
        case 11
            drmS=findrespmat(r0,indBPM,indSkewCor,kval,'PolynomA',1,2,'finddispersion6Err',indrfc,alpha,delta,inCOD);
            drmS{1}=drmS{1}./kval;
            drmS{2}=drmS{2}./kval;
            drmS{3}=drmS{3}./kval;
            drmS{4}=drmS{4}./kval;
            disp(' --- computed H dispersion RM skew quad from model --- ');
            
            % store data
            ModelRM.DispSCor=drmS;
            
        case 12
            
            drmT=findrespmat(r0,1,indQuadCor,kval,'PolynomB',1,2,'gettunechromatlinopt',inCOD);
            drmT{1}=drmT{1}./kval;
            drmT{2}=drmT{2}./kval;
            drmT{3}=drmT{3}./kval;
            drmT{4}=drmT{4}./kval;
            
            ModelRM.TuneQCor=drmT;
            
            disp(' --- computed tune and chrom RM normal quad from model --- ');
            
        case 13
            
            ModelRM.misal=1e-5;
            ModelRM.DK=1e-3;
            QM=findrespmatmisal(r0,indBPM,indQuadCor,ModelRM.DK,'PolynomB',1,2,ModelRM.misal,inCOD);
            QM{1}=QM{1}./kval;
            QM{2}=QM{2}./kval;
            QM{3}=QM{3}./kval;
            QM{4}=QM{4}./kval;
            
            ModelRM.BBAQuad=QM;
            
            disp(' --- computed orbit introduced by normal quad DK offseted by 100um --- ');
            
            
        otherwise
    end
end

ModelRM.kval=kval;
ModelRM.delta=delta;

return
