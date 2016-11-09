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


disp(['Going to tune: ' num2str(t0,'%2.5f, ')]);
rerr = atfittune(rerr,WPtune,qfidx,qdidx);
rerr = atfittune(rerr,WPtune,qfidx,qdidx);

[b]=atlinopt(rerr,0,1:length(rerr)+1);
tf=b(end).mu/2/pi;
disp(['Final tune: ' num2str(tf,'%2.5f, ')]);

%if ~isempty(find((tf-t0)>=0.1,1)) || isnan(tf)
%    disp('Corrected tune is >0.1 far from nominal')
    rerr=atmatchtunedelta(rerr,WPtune,{qfidx, qdidx});
%end

%
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
