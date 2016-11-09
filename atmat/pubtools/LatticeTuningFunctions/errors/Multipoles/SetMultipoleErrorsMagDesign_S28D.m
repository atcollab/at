function [....
    ringmultipoles,.... % ring with multipole errors
    bnsyst,....         % normal systematic multipole components applied
    bnrand,....         % normal random multipole components applied
    anrand....          % skew   random multipole components applied
    ]=SetMultipoleErrorsMagDesign_S28D(...
    ring,...                %lattice without multipoles errors
    systematicerrorflag,... %set systematic errors
    randomerrorflag,...      %set random errors
    scalerandmulterr,... % scale random multipole errors
    dipquaderrors)       % set also dipole and quadrupole multipole errors
% function ringmultipoles=SetMultipoleErrorsMagDesign_S28BINJ(ring)
% assigns S28BINJ multipole errors
% [....
%     ringmultipoles,.... % ring with multipole errors
%     bnsyst,....         % normal systematic multipole components applied
%     bnrand,....         % normal random multipole components applied
%     anrand....          % skew   random multipole components applied
%     ]=SetMultipoleErrorsMagDesign_S28BINJ(...
%     ring,...                %lattice without multipoles errors
%     systematicerrorflag,... %set systematic errors
%     randomerrorflag...      %set random errors
%     scalerandmulterr) scale random multipole errors 
%           20/30 for bulk
%           40/30 for laminated
%
% update: 20/02/2015 (missing DL (rand), OF, DQ1, OF rand and syst)
% update: July/2015 (updated multipole values)
% update: September/2015 (updated multipole values, added DQ1 DQ2 random assumed DQ1=DQ2 and anr=bnr)
%
%
%see also: AssignFieldErr TruncatedGaussian


if nargin==1
    disp('Setting All multipoles')
    systematicerrorflag=1;
    randomerrorflag=1;
end

warning('Injection magnets as cell magnets!');

magfamsnames={...
    'QF1',...
    'QD2',...
    'DL1_1',...
    'DL1_2_5',...
    'QD3',...
    'SD1',...
    'QF4',...
    'SF ',...
    'OF ',...
    'QD5',...
    'QF6',...
    'DQ1',...
    'QF8',...
    'DQ2'};

magindexes={... % magnet indexes, in order of apperance in the lattice cell
    findcells(ring,'FamName','QF1\w*'),...% [A-E] is to exclude the Injection magnets
    findcells(ring,'FamName','QD2\w*'),...
    findcells(ring,'FamName','DL\w*_1'),...
    findcells(ring,'FamName','DL\w*_[2-5]'),...
    findcells(ring,'FamName','QD3\w*'),...
    findcells(ring,'FamName','S[DJ]1\w*'),...
    findcells(ring,'FamName','QF4[A-Z]'),...
    findcells(ring,'FamName','S[FIJ]2\w*'),...
    findcells(ring,'FamName','O[FIJ]\w*'),...
    findcells(ring,'FamName','QD5\w*'),...
    findcells(ring,'FamName','QF6\w*'),...
    findcells(ring,'FamName','DQ1\w*'),...
    findcells(ring,'FamName','QF8\w*'),...
    findcells(ring,'FamName','DQ2\w*'),...
    };

onoff=ones(1,length(magindexes));% all on by default

onoff=onoff*1;%*1e-2; % scale all down

if nargin==4
% flag to switch on or off dipole and quadrupole multipole errors.
dipquaderrors=0;
end

offmult=find(onoff==0);
if ~isempty(offmult)
    warning(['Setting to zero SOME multipoles!']);
    for iii=1:length(offmult)
    disp(['setting to zero ' magfamsnames{offmult(iii)} ' multipoles!']);
    end
end
disp(['Scaling multipoles by '  num2str(onoff(1))]);
    
%% assign  known field errors.
%%%%% SYSTEMATIC

% qf1
qf1bn=[  
    0
10000
0
0
0
1.10934
0
0
0
-5.18658
0
0
0
4.45567
0
0
0
-0.922755
0
0
0
-0.0349446
0
0
0
0.0401289
0
0
0
-0.00454002
0
0
   ];

% qd2 qf4 qd5 
qd2qf4qd5bn=[
   0
10000
0
0
0
1.05689
0
0
0
-5.0583
0
0
0
4.30472
0
0
0
-0.900747
0
0
0
-0.032636
0
0
0
0.0392085
0
0
0
-0.00455189
0
0
 ];

% qd3
qd3bn=[
    0
10000
0
0
0
0.213039
0
0
0
-5.38733
0
0
0
4.24571
0
0
0
-0.869278
0
0
0
-0.0369974
0
0
0
0.0389081
0
0
0
-0.00427762
0
0
];

% qf6

qf6bn=[
    0
10000
0
0
0
-0.141575
0
0
0
0.37436
0
0
0
-0.416586
0
0
0
0.0331539
0
0
0
-0.00202149
0
0
0
0.000173277
0
0
0
-0.000015608
0
0
];

% qf8 
qf8bn=[
0
10000
0
0
0
-0.0286179
0
0
0
0.507304
0
0
0
-0.518238
0
0
0
0.0519398
0
0
0
-0.00397534
0
0
0
0.000326059
0
0
0
-2.60953E-05
0
0    
];



dq1bn=[%%%% relative to DIPOLE!
    10000
4499.22*dipquaderrors*0%% considered elesewhere
1.3194
-3.52302
-1.37691
9.23649
1.06358
-10.1486
-3.14353
2.6593
2.20996
-0.0810461
-1.00213
-0.534322
0.0956301
0.230463
0.0973052
-0.028557
-0.0408101
-0.0100374
0.0117989
0.016608
-0.0180537
-0.0091815
-0.0033418
0.00478292
0.017508
0.0104226
0.0563036
0.0103595
0.140688
0.0239935
];

dq2bn=[ % alfresco 09/07/2015
10000
5532.12*dipquaderrors*0%% considered elesewhere
-0.638178
25.8143
5.20062
-17.2604
-5.60715
-3.13825
0.700903
0.594722
-0.0168115
-0.0909287
-0.099857
0.126364
0.102901
0.0660892
];

 
dl25bn=[
 10000
-0.114587528*dipquaderrors
-10.90737749
0.000826888
-9.553706569
0.001844252
40.27705914
0.000985885
-39.06399784
-0.02128641
12.18333403
0.042352407
4.097949452
-0.048263899
-8.461398324
0.033963307
7.417928542
0.006952543
-5.297735568
-0.065757173
3.677007672
0.1190935
-2.698133257
-0.170708161
2.052958641
0.210328708
-1.477094684
-0.236478376
0.896419778
0.241228991
-0.363207639
-0.212722702
-0.048637051
0.171068182
0.277425745
-0.129839868
-0.306962067
0.124492321
0.170936777
 ];

dl1bn=[
    10000
-0.120136117*dipquaderrors
-18.78769519
-0.04104429
-9.229796041
0.175741454
39.31519315
-0.047200252
-37.91431184
0.057877333
11.99815634
0.029617722
3.88913231
-0.116902012
-8.01975352
0.115098705
7.090914685
-0.027823358
-5.007313249
-0.043452941
3.516510125
0.127367135
-2.520179322
-0.131847244
1.922521376
0.081749897
-1.303745402
-0.008143467
0.728894405
-0.034517217
-0.149670856
0.071987846
-0.263699771
-0.148112903
0.553683203
0.162973312
-0.621075716
-0.14177974
0.525840645
];

 

sfdbn=[
    0,...SD1
    0,...SD1
    10000,...SD1
    0,...SD1
    0,...SD1
    0,...SD1
    0,...SD1
    19.996,...SD1
    0,...SD1
    0,...SD1
    0,...SD1
    0,...SD1
    0,...SD1
    0,...SD1
    -7.5048,...SD1
    0,...SD1
    0,...SD1
    0,...SD1
    0,...SD1
    0,...SD1
    -0.0510856,...SD1
    ]';

ofbn=[
    0
    0
    0
    10000
    0
    0
    0
    0
    0
    0];




%% build cellarray of multipole components

bnsyst={...% multipole field errors, in order of apperance in the lattice cell
    [...QF1
    qf1bn]*1e-4*onoff(1),...QF1
    [...QD2
    qd2qf4qd5bn]*1e-4*onoff(2),...QD2
    [...DL1_1
    dl1bn]*1e-4*onoff(3),...DL1
    [...DL1_2-5
    dl25bn]*1e-4*onoff(4),...DL1
    [...QD3
    qd3bn]*1e-4*onoff(5),...QD3
    [...SD1
    sfdbn]*1e-4*onoff(6),...SD1
    [...QF4
    qd2qf4qd5bn]*1e-4*onoff(7),...QF4
    [...SF2
    sfdbn]*1e-4*onoff(8),...SF2
    [...OF1
    ofbn]*1e-4*onoff(9),...OF1
    [...QD5
    qd2qf4qd5bn]*1e-4*onoff(10),...QD5
    [...QF6
    qf6bn]*1e-4*onoff(11),...QF6
    [...DQ1 
    dq1bn]*1e-4*onoff(12),...DQ1 
    [...QF8
    qf8bn]*1e-4*onoff(13),...QF8
    [...DQ2
    dq2bn]*1e-4*onoff(14),...DQ2 
    };


%% RANDOM
% 30 um

% normal
qf1qd2qd3qf4qd5_bnr=[
    0
    0
    3.072975926
2.54644786
2.219625361
1.830077863
1.525070225
1.224277182
1.037263493
0.864779645
0.723679911
0.54251811
0.473629212
0.372401109
]*scalerandmulterr;
qf1qd2qd3qf4qd5_bnr(2)=1e4; % ref multipole

qf1qd2qd3qf4qd5_anr=[
    0
    0
    3.072975926
2.54644786
2.219625361
1.830077863
1.525070225
1.224277182
1.037263493
0.864779645
0.723679911
0.54251811
0.473629212
0.372401109
    ]*scalerandmulterr;

dq1bnr=[
    0
    10000
    9.360035
    6.231414
    4.040027
    2.390795
    1.598305];

dq2bnr=[
    0
    10000
    9.360035
    6.231414
    4.040027
    2.390795
    1.598305];

dlbnr=[
    10000
    0
    0
    0
    0
    0];

sextbnr=[
    0.0022*dipquaderrors
0.0017*dipquaderrors
0
7.2330e-04
3.9599e-04
1.4986e-04
4.2453e-05
1.2127e-05
2.7490e-05
3.5432e-05
3.1562e-05
2.0812e-05
1.5015e-05
8.5440e-06
3.9749e-06
2.8271e-06
1.3636e-06
4.1152e-07
]*1e4*30/50*scalerandmulterr;% computed for 50um
sextbnr(3)=1e4;

octbnr=[
    0
    0
    0
    10000
    0
    0
    0
    0
    0];

% skew

qf6qf8_bnr=[
 0   
0
4.803458371
1.910276957
1.055734675
0.588073151
0.312742308
0.175288289
0.101114708
0.064747269
0.038242921
0.022608539
0.01461673
0.008751506
]*scalerandmulterr;
qf6qf8_bnr(2)=1e4;

qf6qf8_anr=[
 0  
0
4.803458371
1.910276957
1.055734675
0.588073151
0.312742308
0.175288289
0.101114708
0.064747269
0.038242921
0.022608539
0.01461673
0.008751506
]*scalerandmulterr;

dq1anr=[
    0
    0
    9.360035
    6.231414
    4.040027
    2.390795
    1.598305];

dq2anr=[
    0
    0
    9.360035
    6.231414
    4.040027
    2.390795
    1.598305];


dlanr=[0
    0
    0
    0
    0
    0];

sextanr=[
    0.0020*dipquaderrors
0.0018*dipquaderrors
0.0013 
7.4844e-04
3.2321e-04
1.6012e-04
3.8491e-05
1.2814e-05
3.4294e-05
3.6640e-05
2.5829e-05
2.2166e-05
1.3582e-05
9.0634e-06
5.5419e-06
2.9239e-06
1.1148e-06
4.3986e-07
]*1e4*30/50*scalerandmulterr;% computed for 50um

octanr=[
    0
    0
    0
    0
    0
    0
    0
    0
    0];


bnrand={...% multipole field errors, in order of apperance in the lattice cell
    [...QF1
    qf1qd2qd3qf4qd5_bnr]*1e-4*onoff(1),...QF1
    [...QD2
    qf1qd2qd3qf4qd5_bnr]*1e-4*onoff(2),...QD2
    [...DL1
    dlbnr]*1e-4*onoff(3),...DL1
    [...DL1_2-5
    dl25bn]*1e-4*onoff(4),...DL1
    [...QD3
    qf1qd2qd3qf4qd5_bnr]*1e-4*onoff(5),...QD3
    [...SD1
    sextbnr]*1e-4*onoff(6),...SD1
    [...QF4
    qf1qd2qd3qf4qd5_bnr]*1e-4*onoff(7),...QF4
    [...SF2
    sextbnr]*1e-4*onoff(8),...SF2
    [...OF1
    octbnr]*1e-4*onoff(9),...OF1
    [...QD5
    qf1qd2qd3qf4qd5_bnr]*1e-4*onoff(10),...QD5
    [...QF6
    qf6qf8_bnr]*1e-4*onoff(11),...QF6
    [...DQ1
    dq1bnr]*1e-4*onoff(12),...DQ1 
    [...QF8
    qf6qf8_bnr]*1e-4*onoff(13),...QF8
    [...DQ2
    dq2bnr]*1e-4*onoff(14),...DQ2 
    };


% skew random multipoles
anrand={...% multipole field errors, in order of apperance in the lattice cell
    [...QF1
    qf1qd2qd3qf4qd5_anr]*1e-4*onoff(1),...QF1
    [...QD2
    qf1qd2qd3qf4qd5_anr]*1e-4*onoff(2),...QD2
    [...DL1
    dlanr]*1e-4*onoff(3),...DL1
    [...DL2-5
    dlanr]*1e-4*onoff(4),...DL2-5
    [...QD3
    qf1qd2qd3qf4qd5_anr]*1e-4*onoff(5),...QD3
    [...SD1
    sextanr]*1e-4*onoff(6),...SD1
    [...QF4
    qf1qd2qd3qf4qd5_anr]*1e-4*onoff(7),...QF4
    [...SF2
    sextanr]*1e-4*onoff(8),...SF2
    [...OF1
    octanr]*1e-4*onoff(9),...OF1
    [...QD5
    qf1qd2qd3qf4qd5_anr]*1e-4*onoff(10),...QD5
    [...QF6
    qf6qf8_anr]*1e-4*onoff(11),...QF6
    [...DQ1
    dq1anr]*1e-4*onoff(12),...DQ1 
    [...QF8
    qf6qf8_anr]*1e-4*onoff(13),...QF8
    [...DQ2
    dq2anr]*1e-4*onoff(14),...DQ2 
    };

%% assign all systematic errors 


if systematicerrorflag
    
    for imag=1:length(magindexes)
        
        if ~isempty(bnsyst{imag}) %&& find(bnsyst{imag})
            % find reference multipole, set it to zero
            bs=bnsyst{imag};
            refMult=find(bs==1);
          
            if isempty(refMult)
                maxord=ring{magindexes{imag}(1)}.MaxOrder+1;
                %ring{magindexes{imag}(1)}
                [~,maxcoef]=max(abs(ring{magindexes{imag}(1)}.PolynomB));
                refMult=min(maxord,maxcoef);
            end
            
            bs(refMult)=0;
            if onoff(imag)
                disp(['Systematic: ' magfamsnames{imag}...
                    ', refmultipole: ' num2str(refMult)...
                    ', std(berr): ' num2str(std(bs))])
                
                % assign field errors
                rho=getcellstruct(ring,'RefRadius',magindexes{imag});
                rho(isnan(rho))=0.013; % default for magnet without refRadius =0.013
                %             find(isnan(rho),1,'first')
                %             unique(rho(~isnan(rho)))
                
                [ring,pb]=AssignFieldErr(ring,magindexes{imag},refMult,rho,bs');%,0*bs
          
                 
               if find(isnan(pb))
                
                 error(['NaN in PolynomB for ' magfamsnames{imag}])
                end
            end
            
         else
            warning(['Missing Systematic multipole data for: ' magfamsnames{imag} ])
         end
    end

end

%% assign all random errors 
if randomerrorflag
 %error('STILL DO NOT KNOW WHAT [arbitrary units] Means in the P folder data!')
    for imag=1:length(magindexes)
        
        if ~isempty(bnrand{imag}) || ~isempty(anrand{imag})
            % find reference multipole, set it to zero
            bs=bnrand{imag};
            refMult=find(bs==1);
            %bsN=bnsyst{imag};
            %refMult=find(bsN==1);
            if isempty(refMult)
                maxord=ring{magindexes{imag}(1)}.MaxOrder+1;
                [~,maxcoef]=max(abs(ring{magindexes{imag}(1)}.PolynomB));
                refMult=min(maxord,maxcoef);
            end
            
            bs(refMult)=0;
            as=anrand{imag};
            as(refMult)=0;
            
            if onoff(imag)
                disp(['Random: ' magfamsnames{imag} ', refmultipole: ' num2str(refMult)])
                indh=magindexes{imag};
                rho=getcellstruct(ring,'RefRadius',magindexes{imag});
                rho(isnan(rho))=0.013; % default for magnet without refRadius =0.013
                Nsig=2;
                
                if find(bs) % only if some exist.
                    ipolb=find(bs);
                    for indpolb=1:length(ipolb)%find(bs); % loop berr to get random errors.
                        bsr(ipolb(indpolb),:)=TruncatedGaussian(bs(ipolb(indpolb)),abs(Nsig*bs(ipolb(indpolb))),length(indh));
                    end
                    
                    if find(as) % only if some exist.
                        ipolb=find(as);
                        for indpolb=1:length(ipolb)%ipolb=find(as); % loop berr to get random errors.
                            asr(ipolb(indpolb),:)=TruncatedGaussian(as(ipolb(indpolb)),abs(Nsig*as(ipolb(indpolb))),length(indh));
                         end
                    else
                        asr=zeros(size(bsr));
                    end
                    
                    for jpolb=1:length(magindexes{imag})  % assign random multipoles to each magnet
                        ring=AssignFieldErr(ring,indh(jpolb),refMult,rho(jpolb),bsr(:,jpolb)',asr(:,jpolb)');
                    end
                end
                clear bsr
                clear asr
            end
        else
            warning(['Missing Random multipole data for: ' magfamsnames{imag} ])
        end
        
    end

end

ringmultipoles=ring;
