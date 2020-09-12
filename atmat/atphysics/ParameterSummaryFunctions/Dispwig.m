function [Ds,Dsp,x] = Dispwig (rho,Lw,Dis0)
% Using Finite Difference method to compute dispersion through Wiggler
% Ds : eta
% Dsp: eta';
% A.Mash'al, Iranian Light Source Facility, 2020-09-03
D0=Dis0;
N=101;
h=Lw/(N-1);
x=0:h:Lw;
f=@(x) rho(x);
g=@(x) f(x).*f(x);
D_old=D0*ones(N,1);
D_old(2:N-1)=D0;
D_new=D_old;
cond=1;
count=0;
while cond==1
    
for i=2:N-1
    D_new(i)=(f(x(i))-((D_old(i+1)+D_new(i-1))/(h^2)))/(g(x(i))-(2/(h^2)));
end

    Err=abs(norm(D_new)-norm(D_old));
    count=count+1;
    
    if Err<1e-9 || count>10000
        cond=0;
    end
    D_old=D_new;
end
fprintf('Iteration = %d , Error= %e\n',count,Err)
Ds=D_new;
Dsp=zeros(N,1);
for i=2:length(x)
Dsp(i)=(Ds(i)-Ds(i-1))/(x(2)-x(1));
end
Dsp(1)=Dsp(end);
for i=2:length(x)
Dsp(i)=(Ds(i)-Ds(i-1))/(x(2)-x(1));
end
end