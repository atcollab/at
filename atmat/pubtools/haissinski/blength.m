function sigz=blength(Q)

sigz = (2/3)^(1/3)/(9*Q + sqrt(3)*sqrt(-4+27*Q^2))^(1/3)...
    + (9*Q + sqrt(3)*sqrt(-4+27*Q^2))^(1/3)/(2^(1/3)*3^(2/3));