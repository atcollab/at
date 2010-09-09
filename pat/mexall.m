
% mex  DriftPass.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

apassfile = dir('*Pass.c');
for ii=1:length(apassfile)
    if strcmp('BndMPoleSymplectic4E2Pass.c',apassfile(ii).name)
        continue
    elseif strcmp('BndMPoleSymplectic4E2RadPass.c',apassfile(ii).name)
        continue
    else
        eval(['mex ' apassfile(ii).name '   CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" ']);
    end
    
end
mex -c atlalib.c
mex -c atphyslib.c

mex  BndMPoleSymplectic4E2Pass.c  atlalib.obj atphyslib.obj CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex  BndMPoleSymplectic4E2RadPass.c  atlalib.obj atphyslib.obj CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
% AperturePass.c                  BndMPoleSymplectic4RadPass.c    IdentityPass.c                  StrMPoleSymplectic4RadPass.c
% BendLinearPass.c                CavityPass.c                    Matrix66Pass.c                  ThinMPolePass.c
% BndMPoleSymplectic4E2Pass.c     CorrectorPass.c                 QuadLinearPass.c
% BndMPoleSymplectic4E2RadPass.c  DriftPass.c                     SolenoidLinearPass.c
% BndMPoleSymplectic4Pass.c       EAperturePass.c                 StrMPoleSymplectic4Pass.c
