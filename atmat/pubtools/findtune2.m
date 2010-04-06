function [tune,numdc]=findtune2(M)
numdc=0;
p=size(M,2);
n=size(M,1);
nu=[0:1:n-1];
han=(0.53836-0.46164*cos(2*pi*nu./(n-1)))';
hanext=han(:,ones(1,p));
 M=M.*hanext;   


res=abs(fft(M));
res(1,:)=0;
[bid,ind]=max(res);


for k=1:size(res,2)
    serie=res(:,k);
    
    if ind(k)==1
    
        tune(k)=0;
        numdc=numdc+1;
    elseif ind(k)==n
    tune(k)=1-1/n;
    
    else 
    a=ind(k);
    interpol=(serie(a)*a+serie(a-1)*(a-1)+serie(a+1)*(a+1))/(serie(a)+serie(a-1)+serie(a+1));
    tune(k)=(interpol-1)/n;
    end
end
    