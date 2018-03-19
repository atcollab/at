function rerr=fittunedelta2fam(rerr,r0)
% rerr=fittunedelta2fam(rerr,r0)
%
% matches the tune of rerr to that of r0.
%
% the 2 quadrupoles families used to correct the tune are marked by field:
% qfidx=findcells(rerr,'ForTuneF');
% qdidx=findcells(rerr,'ForTuneD');
%
% if atfittune fails a second attempt is made using atmatchtunedelta
%
%see also: atfittune atmatchtunedelta

disp('Tune Matching')
[b]=atlinopt(r0,0,1:length(r0)+1);
t0=b(end).mu/2/pi;
disp(['Nominal tune: ' num2str(t0,'%2.5f, ')]);
WPtune=t0;

qfidx=findcells(rerr,'ForTuneF');
qdidx=findcells(rerr,'ForTuneD');

rerr0=rerr;% inital tune lattice

[b]=atlinopt(rerr,0,1:length(rerr)+1);
te=b(end).mu/2/pi;
disp(['Initial tune: ' num2str(te,'%2.5f, ')]);

% % match also integer part of the tune
% disp(['Going to tune: ' num2str(t0+(te-WPtune)*0.5,'%2.5f, ')]);
%
% rerr=atmatchtunedelta(rerr,t0+(te-WPtune)*0.5,{qfidx, qdidx});
%
% [b]=atlinopt(rerr,0,1:length(rerr)+1);
% ti=b(end).mu/2/pi;
% disp(['Intermediate tune: ' num2str(ti,'%2.5f, ')]);
%
% rerr=atmatchtunedelta(rerr,WPtune,{qfidx, qdidx});
%

%% loop to find correctible tunes (non integer)
Dtunetest=[0,0];
Dtune = 0.05;
modk=atgetfieldvalues(r0,[qdidx qfidx],'PolynomB',{1,2});
errk=atgetfieldvalues(rerr,[qdidx qfidx],'PolynomB',{1,2});
rerrt= rerr; 
it = +4;

frac=@(x)(x-floor(x));

tef=frac(te);

% impose change if tune close to integer or half integer
if tef(1)<0.1
    disp(' Qh<0.1')
    te(1) = NaN;
end
if tef(2)<0.1
    disp(' Qv<0.1')
    te(2) = NaN;
end
if tef(1)>0.9
    disp(' Qh > 0.9')
    te(1) = NaN;
end
if tef(2)>0.9
    disp(' Qv > 0.9')
    te(2) = NaN;
end
if tef(1)>0.45 && tef(1)<0.55
    disp(' 0.45<Qh<0.55')
    te(1) = NaN;
end
if tef(2)>0.45 && tef(2)<0.55
    disp(' 0.45<Qv<0.55')
    te(2) = NaN;
end

while  ~isempty(find(isnan(te)==1,1)) && it>-4
    
    Dtunetest(isnan(te))=Dtune*it;
    it = it-1;
    disp(['Initial tune is NaN, by ' num2str(Dtunetest) ', to find a readable tune'])
    tt = WPtune + Dtunetest;
    r0p1 = atfittune(r0,tt-floor(tt),qfidx,qdidx);
    r0p1 = atfittune(r0p1,tt-floor(tt),qfidx,qdidx);
    dtuk=atgetfieldvalues(r0p1,[qdidx qfidx],'PolynomB',{1,2});
    
    % set Delta tune of 0.1 in both planes, to search for readable tune
    % value (move lattice with errors from integer, half integer resonances)
    rerrt = atsetfieldvalues(rerr,[qdidx qfidx],'PolynomB',{1,2},errk+(dtuk-modk));
    
    % try tune with errors again
    [b]=atlinopt(rerrt,0,1:length(rerrt)+1);
    te=b(end).mu/2/pi;
    disp(['Initial tune: ' num2str(te,'%2.5f, ')]);
end

rerr = rerrt;

disp(['Going to tune: ' num2str(t0,'%2.5f, ')]);
if abs(te-WPtune) > [0.1 0.1] % slowly go to nominal tune
    
    for fracval=(0.75:-0.25:0)
        tt = WPtune+(te-WPtune)*fracval;
        
        rerrt = atfittune(rerrt,frac(tt),qfidx,qdidx);
        
        % if tune set ok store improved lattice
        [b]=atlinopt(rerrt,0,1:length(rerrt)+1);
        te=b(end).mu/2/pi;
        
        disp(['Intermediate tune: ' num2str(te,'%2.5f, ')]);
        
        if isempty(find(isnan(te)==1,1))
            rerr = rerrt;
        end
    end
    
else
    % go in one shot
    
    rerrt = atfittune(rerrt,WPtune,qfidx,qdidx);
    
    disp(['Intermediate tune: ' num2str(te,'%2.5f, ')]);
    
    % if tune set ok store improved lattice
    [b]=atlinopt(rerrt,0,1:length(rerrt)+1);
    te=b(end).mu/2/pi;
    if isempty(find(isnan(te)==1,1))
        rerr = rerrt;
    end
end

[b]=atlinopt(rerr,0,1:length(rerr)+1);
tf=b(end).mu/2/pi;
disp(['Final tune: ' num2str(tf,'%2.5f, ')]);

% grant that integer tune is correct
%disp(['check integer tune ok, match full tunes']);
%rerr=atmatchtunedelta(rerr,WPtune,{qfidx, qdidx});

%%
if ~isempty(find(isnan(tf)==1,1))
    disp('Corrected tune is NaN')
    rerr=rerr0;
    error('Corrected tune is NaN')
end

if ~isempty(find(imag(tf)~=0,1))
    disp('Corrected tune is not real')
    rerr=rerr0;
    error('Corrected tune is not real')
end

if ~isempty(find((tf-t0)>=0.1,1))
    disp('Corrected tune is >0.1 far from nominal')
    rerr=rerr0;
    error('Corrected tune is >0.1 far from nominal')
end

end
