function [h3,h4,h3chro,adts]=RDTfluctuationIndicator(ring,varargin)
%RDTFLUCTUATIONINDICATOR quantitative representation of RDT fluctuations
%   This function calls computeRDTfluctuation(ring, varargin) to compute
%   RDT fluctuations, and provides one example to quantitatively represents 
%   the RDT fluctuations.
%
% [h3,h4,h3chro,adts]=RDTfluctuationIndicator(ring,varargin)
%
%   ring is the AT lattice
%   The additional argument:
%     nslices: number of slices of each sextupole, which affects the
%              crossing terms. default: 4.
%     method:  method used to aggregate the RDT fluctuations into a scalar
%              indicator. 'rms' (default) uses the root-mean-square of the
%              RDTs; 'avg' uses their mean (the original algorithm).
%     dp:      momentum deviation. default: 0.0.
% 
%   h3 and h4 quantitatively represents 3rd- and 4th-order geometric 
%     RDT fluctuations, respectively. 
%
%   h3 + w * h4 can be used to represents geometric RDT fluctuations.
%     The coefficient w can be estimated by action of particle at the
%     anticipated DA, w ~ (2J_x)^0.5 = x / betax^0.5,  usually 0.01 is OK.
% 
%   h3chro is the fluctuation of 3rd-order chromatic RDTs, 
%     defined similarly to h3. 
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
%   [2] B. Wei, Z. Bai, G. Feng, P. Brunelle, A. Loulergue, and L. S. Nadolski, Phys. Rev. Accel. Beams 29, 083404 (2026)
%
[dp,varargin]=getoption(varargin,'dp',0.0);
[method,varargin]=getoption(varargin,'method','rms');

% Select the aggregation function according to 'method'
% ('rms' is implemented as sqrt(mean(x.^2)) to avoid a toolbox dependency)
switch lower(method)
    case 'avg'
        aggregate = @(x) mean(x);
    case 'rms'
        aggregate = @(x) sqrt(mean(x.^2));
    otherwise
        error('RDTfluctuationIndicator:method', ...
              'method must be ''avg'' or ''rms'', got ''%s''.', method);
end

[RDT, ~, natural] = computeRDTfluctuation(ring, 'nperiods', 1, 'dp', dp, varargin{:});
h3 = sqrt(aggregate(natural.f21000)^2 + aggregate(natural.f30000)^2 ...
           + aggregate(natural.f10110)^2 + aggregate(natural.f10020)^2 ...
           + aggregate(natural.f10200)^2);
if dp ~= 0
    % remove misleading information. 
    % the off-momentum RDTs directly describe the nonlinear dynamics for
    % the given dp case. So we don't use the chromatic terms which
    % describe the deviation from the on-momentum case.
    h3chro = NaN;
else
    h3chro = sqrt(aggregate(natural.f20001)^2 + aggregate(natural.f00201)^2 ...
           + aggregate(natural.f10002)^2);
end
h4 = sqrt(aggregate(natural.f31000)^2 + aggregate(natural.f40000)^2 ...
           + aggregate(natural.f20110)^2 + aggregate(natural.f11200)^2 ...
           + aggregate(natural.f20020)^2 + aggregate(natural.f20200)^2 ...
           + aggregate(natural.f00310)^2 + aggregate(natural.f00400)^2);
adts = sqrt(RDT.dnux_dJx^2 + RDT.dnux_dJy^2 + RDT.dnuy_dJy^2);
end