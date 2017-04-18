function z = atskewquad(fname,L , Qs , method)
%atskewquad(Fname, L, Qs, method)
    z = atmultipole(fname,L, [0 Qs 0 0], [0 0 0 0], method);