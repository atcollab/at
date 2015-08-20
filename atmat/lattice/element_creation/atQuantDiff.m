function elem=atQuantDiff(fname,arg)
%atQuantDiff creates a quantum diffusion element
%
%ELEM=ATQUANTDIFF(FAMNAME,DIFFMAT) uses the given diffusion matrix
%   FAMNAME:   family name
%   DIFFMAT:   Diffusion matrix
%
%ELEM=ATQUANTDIFF(FAMNANE,RING) computes the diffusion matrix of the ring
%   FAMNAME:   family name
%   RING:      lattice without radiation
%
%See also quantumDiff

if iscell(arg)
    [ring2,radindex]=atradon(arg);
    dmat=quantumDiff(ring2,radindex);
else
    dmat=arg;
end

elem=atbaselem(fname, 'QuantDiffPass', 'Class', 'QuantDiff', 'Lmatp' , lmatp(dmat));

    function lmatp = lmatp(dmat)
        %lmat does Cholesky decomp of dmat unless diffusion is 0 in
        %vertical.  Then do chol on 4x4 hor-long matrix and put 0's
        %in vertical
        try
            lmat66 = chol(dmat);
        catch
            lm=[chol(dmat([1 2 5 6],[1 2 5 6])) zeros(4,2);zeros(2,6)];
            lmat66=lm([1 2 5 6 3 4],[1 2 5 6 3 4]);
        end
        lmatp=lmat66';
    end
end


