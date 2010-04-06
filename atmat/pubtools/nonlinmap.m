function [nuh,nuv,X,Y,orh,orv,NU,KSI]=nonlinmap(machine,X,Y,nbturn,type,DP)
%Computes the frequency map on the bases of a FFT analysis on tracked
%orbits using the AT routines for traking.
%INPUTS: MACHINE is the at structure representing the accelerator with no
%time dependent fields (no cavities
%X is a vector containing the initial conditions for the horizontal
%position
%Y is a vector containing the initial position for the vertical position
%the number of initial considions tested is size(X)*size(Y). for each X all
%different Y will be tested (ans vis-versa)
%(if type=1) or the momentum deviation (if type=2).
%NBTURN is the number of turns for the traking.Traking takes place for
%2*NBTURN and each half is processed separately.
%TYPE should be 1 for an XY map, 2 for a XDP map.
%DP relevant only if TYPE=1. Defines the momentum deviation for which the
%XY map is calculated.

%OUTPUTS:NUH and NUV are matrices of sizes: size(x)*size(Y) 4
% each lines corresponds to a given initial conditions.
% The initial conditions are aligned in the following way:(X1,Y1)
% (X2,Y1)...(Xn,Y1) (X1,Y2) (X2,Y2)...(Xn,Y2).....(X1,Yn)...(Xn,Yn).
%the first row is the tune computed on the first NBTURNS turns.
%the second row is te tune for the next NBTURN turns.
%the third row is the square of the difference of the first two rows
%the fourth row is the DC orbit.
%H refers to horizontal, V to vertical.
% The initial conditions are aligned in the following way:(X1,Y1)
% (X2,Y1)...(Xn,Y1) (X1,Y2) (X2,Y2)...(Xn,Y2).....(X1,Yn)...(Xn,Yn).


mode=1;
nbturn=nbturn*2;
machine=atreduce(machine);
[LinData,NU, KSI] = LINOPT(machine,0,1);
tune=NU
chroma=KSI
chromanormaesrf=KSI./(NU+[36.44 13.39])
resp=0;
nbpart=numel(X)*numel(Y)  

%creates the bunch of initial particuls. Particules are having a given
%angle in order to avoid null orbits
Rin=zeros(6,nbpart);
Rin([2 4],:)=1E-6*0.32E-17;
[y,x]=meshgrid(Y,X);

%Fills in the initial conditions in the particule bunch 
if type==1
    Rin(5,:)=DP;
    Rin([1 3],:)=[x(:)';y(:)'];
    end 
  
    if type==2
    disp=findorbit4(machine,0.01)
    Rin([1 5],:)=[x(:)';y(:)']; 
    
    
    end
%puts a first test particule with perferct initial condition to avoit an AT
%behaviour that kills the bunch if the first particule is lost
    Rintest=zeros(6,1);
Rin=[Rintest Rin];        
nbpart=nbpart+1;        

%calculates and adds the chromatic orbit to the initial conditions
disp=findorbit4(machine,0.01);    
Rchroma=(Rin(5,:)'*disp([1 3])'*100)';
Rin([1 3],:)=Rin([1 3],:)+Rchroma;
rc=Rchroma(:);
rc2=rc(:,ones(1,nbturn));   
Routcor=(reshape(rc2(:),2,[]));

%tracking routine
for k=1:nbturn         
    [out,loss]=ringpass(machine,Rin,1);
    out(6,:)=0;
    Rfin(:,(k-1)*nbpart+1:(k)*nbpart)=out;
    Rin=out;
end

%substracts the chromatic orbit calculated 
Rfin([1 3],:)=Rfin([1 3],:)-Routcor;
%extract the horizontal and vertical position and momentum from the output matrix
    orh=(reshape(Rfin(1,:),nbpart,nbturn))';
    ph= (reshape(Rfin(2,:),nbpart,nbturn))';
    orv=(reshape(Rfin(3,:),nbpart,nbturn))';
    pv=(reshape(Rfin(4,:),nbpart,nbturn))';
    %kills the usless first particule
    orh(:,1)=[];
    orv(:,1)=[];
ph(:,1)=[];
pv(:,1)=[];    
    
loss(1)=[];


%separates the orbits of the non lost particules in two halfs corresponding
%to the first ans second NBTURNS 
notloss=not(loss);
nbpart=nbpart-1;
orhred=reshape(orh(:,notloss),nbturn/2,[]);
orvred=reshape(orv(:,notloss),nbturn/2,[]);
phred=reshape(ph(:,notloss),nbturn/2,[]);
pvred=reshape(pv(:,notloss),nbturn/2,[]);
%calculates the mean orbit 
meanh=(mean(orh(:,notloss),1))';
meanv=(mean(orv(:,notloss),1))';
maxh=(max(abs(orh(:,notloss)),[],1))';
maxv=(max(abs(orv(:,notloss)),[],1))';
%finds the tunes in both planes for each orbit
[resh,numdch]=findtune2(orhred-i.*phred);
[resv,numdcv]=findtune2(orvred-i.*pvred);
%Fills in the matrices for Nuh and Nuv
nuh(1:nbpart,1:5)=NaN;
nuv(1:nbpart,1:5)=NaN;
deltanu(1:nbpart)=NaN;


tuneh=resh';
tunev=resv';

nuh(notloss,1)=(tuneh(1:2:numel(tuneh)));
nuh(notloss,2)=(tuneh(2:2:numel(tuneh)));
nuh(notloss,3)=(nuh(notloss,1)-nuh(notloss,2)).^2;
nuv(notloss,1)=(tunev(1:2:numel(tunev)));
nuv(notloss,2)=(tunev(2:2:numel(tunev)));
nuv(notloss,3)=(nuv(notloss,1)-nuv(notloss,2)).^2;
nuh(notloss,4)=meanh;
nuv(notloss,4)=meanv;
nuh(notloss,5)=maxh;
nuv(notloss,5)=maxv;
lostpart=sum(loss)