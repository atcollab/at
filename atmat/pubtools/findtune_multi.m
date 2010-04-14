function [tune,spectrum]=findtune_multi(pos,method)

if nargin < 2, method=3; end

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
[vmax,rmax]=max(abs(spectrum(1:nt2,:)));
tune=(rmax-1)/nturns;

case 2
methname='interpolation';
pos2=pos-posm(ones(nturns,1),:);
spectrum=fft(pos2);
[vmax,rmax]=max(abs(spectrum(1:nt2,:)));
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
[vmax,rmax]=max(abs(spectrum(1:nt2,:)));
rmax=rmax+(rmax==1);
kmax=sub2ind([nturns nparts],rmax,1:nparts);
back=(spectrum(kmax-1) > spectrum(kmax+1));
k1=kmax-back;
k2=k1+1;
v1=abs(spectrum(k1));
v2=abs(spectrum(k2));
tune=(rmax-back-1 +((2*v2-v1)./(v1+v2)))/nturns;
%tune2=(rmax-back-1)/nturns + asin(phi(v1,v2,cos(2*pi/nturns))*sin(2*pi/nturns))/2/pi;
%disp(['method 4 tune: ' num2str(mean2(tune2')) ' (rms: ' num2str(std2(tune2')) ')']);

case 4
methname='naff';
pos2=pos-posm(ones(nturns,1),:);
[tune,amplitude,phase]=specfran(pos2);
tune=abs(tune);

end
tune(wrong)=NaN;
errmax=2.5*std2(tune');
keep=(abs(tune-mean2(tune'))<=errmax);
reject=find(~(keep | wrong));
for bpm=reject
    [bname,kdx]=srbpmname(bpm);
    disp(['reject ' bname ' (' num2str(kdx) ')']);
end
disp(sprintf('%20s tune:%g (rms:%g)',methname, mean2(tune(keep)'),std2(tune(keep)')));

function vv=phi(a,b,c)
d1=c*(a+b);
delt=d1.*d1 - 2*a.*b.*(2*c*c-c-1).*a.*a - b.*b - 2*a.*b*c;
vv=(-(a+b*c).*(a-b) + b.*sqrt(delt))./(a.*a + b.*b  +2*a.*b*c);


function w=hann_window(n)
w=0.5*(1-cos(2*pi*(0:n-1)'/n));
