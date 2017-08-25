function z = skewquad(fname,L , Qs , method)
%skewquad(Fname, L, Qs, method)
    z = multipole(fname,L, [0 Qs 0 0], [0 0 0 0], method);