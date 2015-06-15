function [deltap,deltam,Rfin, loss,X0l] = atmomap(ring, Nturn)
%ATMOMAP			find momentum aperture at start of ring
%
%[deltap,deltam] = atmomap(THERING, Nturn)

%dpp0 = [-0.04:0.001:0.04];
%dpp0 = [-0.04:0.005:0.04];
%dpp0 = [-0.04:0.002:-0.02, 0, 0.02:0.002:0.04];
%dpp0 = [-0.05:0.001:-0.02, 0, 0.02:0.001:0.05];
dpp0 = [-0.05:0.05:-0.01, 0, 0.01:0.05:0.05];
%dpp0=0
X0l = zeros(6, length(dpp0));
X0l(5,:) = dpp0;
%tic
%X0l=[0 0 0 0 -0.03 0]';
[Rfin, loss] =ringpass(ring,X0l,Nturn);
%Rfin
%toc
[tmp, indxzero] = find(dpp0==0);
taglossp = min([length(dpp0) indxzero-1+find(loss(indxzero:end))]);
taglossm = max([1 find(loss(1:indxzero))]);

deltap = dpp0(taglossp);
deltam = dpp0(taglossm);
