function J = DampJfromTau(dampingtimes)
%
%   J = DampJfromTau([taux,tauy,tauz])
%
%   J is the vector [Jx, Jy, Jz]
%

alphax=1/dampingtimes(1);
alphay=1/dampingtimes(2);
alphaz=1/dampingtimes(3);

alphatot=alphax+alphay+alphaz;

Jx=4*alphax/alphatot;
Jy=4*alphay/alphatot;
Jz=4*alphaz/alphatot;

J=[Jx, Jy, Jz];
end