function phase = phis(U0MeV,VrfMV)
% this function returns the synchronous phase in radians
% input:
% VrfMV is the RF voltage in MV
% U0MeV is energy loss per turn in MeV

phase=pi - asin(U0MeV./VrfMV);