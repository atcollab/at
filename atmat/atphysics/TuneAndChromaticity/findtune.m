function [tune,spectrum]=findtune(pos,varargin)
%FINDTUNE   get the tune value from turn by turn positions
%
%TUNE=FINDTUNE(POS,METHOD)
%
%POS:       Tune-by-turn particle position
%
%OPTIONAL INPUTS:
%METHOD:    Method for tune determination:
%               1: Highest peak in fft
%               2: Interpolation on fft results
%               3: Windowing + interpolation (default)
%               4: NAFF
%ANG:       Tune-by-turn particle angle (used only for NAFF)
%
%[TUNE,SPECTRUM]=FINDTUNE(...) Also returns the fft


[method,ang]=getargs(varargin,{3,[]});

if method==4 && isempty(ang)
    warning('NAFF method requires particle angles as 3rd input. Switching to method 3');
    method=3;
end


nturns=size(pos,1);
nparts=size(pos,2);
nt2=fix(nturns/2);

posm=mean(pos);
wrong=~isfinite(posm);

switch method
    case 1
        methname='highest peak';
        pos2=pos-posm(ones(nturns,1),:);
        spectrum=fft(pos2);
        [vmax,rmax]=max(abs(spectrum(1:nt2,:))); %#ok<ASGLU>
        tune=(rmax-1)/nturns;
        
    case 2
        methname='interpolation';
        pos2=pos-posm(ones(nturns,1),:);
        spectrum=fft(pos2);
        [vmax,rmax]=max(abs(spectrum(1:nt2,:))); %#ok<ASGLU>
        rmax=rmax+(rmax==1);
        kmax=sub2ind([nturns nparts],rmax,1:nparts);
        back=(spectrum(kmax-1) > spectrum(kmax+1));
        k1=kmax-back;
        k2=k1+1;
        v1=abs(spectrum(k1));
        v2=abs(spectrum(k2));
        tune=(rmax-back-1 +(v2./(v1+v2)))/nturns;
        
    case 3
        methname='window + interp.';
        w=hann_window(nturns);
        pos2=(pos-posm(ones(nturns,1),:)).*w(:,ones(1,nparts));
        spectrum=fft(pos2);
        [vmax,rmax]=max(abs(spectrum(1:nt2,:))); %#ok<ASGLU>
        rmax=rmax+(rmax==1);
        kmax=sub2ind([nturns nparts],rmax,1:nparts);
        back=(spectrum(kmax-1) > spectrum(kmax+1));
        k1=kmax-back;
        k2=k1+1;
        v1=abs(spectrum(k1));
        v2=abs(spectrum(k2));
        tune=(rmax-back-1 +((2*v2-v1)./(v1+v2)))/nturns;
        %tune2=(rmax-back-1)/nturns + asin(phi(v1,v2,cos(2*pi/nturns))*sin(2*pi/nturns))/2/pi;
        %disp(['method 4 tune: ' num2str(mean2(tune2')) ' (rms: ' num2str(std(tune2')) ')']);
    case 4
        methname='naff';
        pos2=pos-posm(ones(nturns,1),:);
        spectrum=nan(nparts,1);
        for iamp=1:nparts
            spectrum(iamp,:)=abs(calcnaff(pos2(:,iamp),ang(:,iamp),1,1))/(2*pi);
        end
        tune=spectrum(:,1)';
end

tune(wrong)=NaN;
errmax=2.5*std(tune');
keep=(abs(tune-mean(tune'))<=errmax);
reject=find(~(keep | wrong));
for bpm=reject
    disp(['reject bpm #' num2str(bpm)]);
end
fprintf('%20s tune:%g (rms:%g)\n',methname, mean(tune(keep)'),std(tune(keep)'));

function vv=phi(a,b,c) %#ok<DEFNU>
d1=c*(a+b);
delt=d1.*d1 - 2*a.*b.*(2*c*c-c-1).*a.*a - b.*b - 2*a.*b*c;
vv=(-(a+b*c).*(a-b) + b.*sqrt(delt))./(a.*a + b.*b  +2*a.*b*c);


function w=hann_window(n)
w=0.5*(1-cos(2*pi*(0:n-1)'/n));
