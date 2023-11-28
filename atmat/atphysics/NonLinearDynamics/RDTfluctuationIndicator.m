function [h3ave,h4ave,h3chroave,adts]=RDTfluctuationIndicator(ring,varargin)
%RDTFLUCTUATIONINDICATOR quantitative representation of RDT fluctuations
%   This function calls computeRDTfluctuation(ring, varargin) to compute
%   RDT fluctuations, and provides one example to quantitatively represents 
%   the RDT fluctuations.
%
% [h3ave,h4ave,h3chroave,adts]=RDTfluctuationIndicator(ring,varargin)
%
%   ring is the AT lattice
%   The additional argument:
%     nslices: number of slices of each sextupole, which affects the
%              crossing terms. default: 4.
% 
%   h3ave and h4ave quantitatively represents 3rd- and 4th-order geometric 
%     RDT fluctuations, respectively. 
%
%   h3ave + w * h4ave can be used to represents geometric RDT fluctuations.
%     The coefficient w can be estimated by action of particle at the
%     anticipated DA, w ~ (2J_x)^0.5 = x / betax^0.5,  usually 0.01 is OK.
% 
%   h3chroave is the fluctuation of 3rd-order chromatic RDTs, 
%     defined similarly to h3ave. 
%
%   adts is the sum of three ADTS terms, which are also used in onlinear optimization.
%   The fluctuations of ADTS terms are not considered.
%   It is calculated here to avoid duplicate computation.
%
% Noting:
%    1.This function provides one example to quantitatively represents the 
%      RDT fluctuations similar to that in Ref.[1]. But in Ref.[1], only
%      the RDTs at the locations of nonlinear magnets are considered. 
%      We think the differences are small and it's more important to 
%      keep the function simple.
%    2.People can call computeRDTfluctuation(ring, varargin) directly,
%      and try other quantitative representations.
%    3.The build-up RDT fluctuation can also be used, see Ref.[1].
%      If the build-up RDT fluctuations are used, it is better to calculate
%      the RDT build-up fluctuations for multiple periods to have
%      better convergence of calculation.
%
% REFERENCE:
%   [1] B. Wei, Z. Bai, J. Tan, L. Wang, and G. Feng, Phys. Rev. Accel. Beams 26, 084001 (2023)
%

[RDT, ~, natural] = computeRDTfluctuation(ring, 'nperiods', 1, varargin);
h3ave = sqrt(mean(natural.f21000)^2 + mean(natural.f30000)^2 ...
           + mean(natural.f10110)^2 + mean(natural.f10020)^2 ...
           + mean(natural.f10200)^2);
h3chroave = sqrt(mean(natural.f20001)^2 + mean(natural.f00201)^2 ...
           + mean(natural.f10002)^2);
h4ave = sqrt(mean(natural.f31000)^2 + mean(natural.f40000)^2 ...
           + mean(natural.f20110)^2 + mean(natural.f11200)^2 ...
           + mean(natural.f20020)^2 + mean(natural.f20200)^2 ...
           + mean(natural.f00310)^2 + mean(natural.f00400)^2);
adts = sqrt(RDT.dnux_dJx^2 + RDT.dnux_dJy^2 + RDT.dnuy_dJy^2);
end