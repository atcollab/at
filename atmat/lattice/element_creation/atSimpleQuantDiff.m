function elem = atSimpleQuantDiff(fname,varargin)
%SimpleQuantDiff creates a simple quantum difusion element
%
%ELEM=SIMPLEQUANTDIFFPASS(FAMNAME,...)
%   FAMNAME:   family name
%
%ELEM=SIMPLEQUANTDIFFPASS(FAMNAME,...,'Betax',BETAX,...)
%   BETAX:   Horizontal beta function. Default: 1.0
%
%ELEM=SIMPLEQUANTDIFFPASS(FAMNAME,...,'Betay',BETAY,...)
%   BETAY:   Vertical beta function. Default: 1.0
%
%ELEM=SIMPLEQUANTDIFFPASS(FAMNAME,...,'Emitx',EMITX,...)
%   EMITX:   Horizontal equilibrium emittance. Default: 0.0
%
%ELEM=SIMPLEQUANTDIFFPASS(FAMNAME,...,'Emity',EMITY,...)
%   EMITY:   Vertical equilibrium emittance. Default: 0.0
%
%ELEM=SIMPLEQUANTDIFFPASS(FAMNAME,...,'Espread',ESPREAD,...)
%   ESPREAD: Equilibrium momentum spread. Default: 0.0
%
%ELEM=SIMPLEQUANTDIFFPASS(FAMNAME,...,'Taux',TAU_X,...)
%   TAU_X: Horizontal damping time. Default: 0.0
%
%ELEM=SIMPLEQUANTDIFFPASS(FAMNAME,...,'Tauy',TAU_Y,...)
%   TAU_Y: Vertical damping time. Default: 0.0
%
%ELEM=SIMPLEQUANTDIFFPASS(FAMNAME,...,'Tauz',TAU_Z,...)
%   TAU_Z: Longitudinal damping time. Default: 0.0
%
%ELEM=SIMPLEQUANTDIFFPASS(FAMNAME,...,'U0',U0,...)
%   U0:     Energy loss [eV]. Default: 0.0

[rsrc,method]=decodeatargs({'SimpleQuantDiffPass'},varargin);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[cl,rsrc]=getoption(rsrc,'Class','SimpleQuantDiff');
[betax,rsrc]      = getoption(rsrc,'Betax',1.0);
[betay,rsrc]      = getoption(rsrc,'Betay',1.0);
[emitx,rsrc]      = getoption(rsrc,'Emitx',0.0);
[emity,rsrc]      = getoption(rsrc,'Emity',0.0);
[espread,rsrc]    = getoption(rsrc,'Espread',0.0);
[taux,rsrc]       = getoption(rsrc,'Taux',0.0);
[tauy,rsrc]       = getoption(rsrc,'Tauy',0.0);
[tauz,rsrc]       = getoption(rsrc,'Tauz',0.0);
[U0,rsrc]         = getoption(rsrc,'U0',0.0);
elem=atbaselem(fname,method,'Class',cl,'Betax',betax,'Betay',betay,...
    'Emitx',emitx,'emity',emity,'Espread',espread,...
    'taux',taux,'tauy',tauy,'tauz',tauz,...
    'U0',U0,rsrc{:});
end