function res = atchapters()

m.id = "atelemcreate";
m.title = "Element creation";
m.list = ["atbaselem", "atelem", "atringparam", "atdrift", "atquadrupole",...
    "atrbend", "atsbend", "atsextupole", "atskewquad", "atmultipole", "atrfcavity",...
    "atmarker", "atmonitor", "ataperture", "atcorrector", "atidtable", "atwiggler",...
    "atdampMatElem", "atsolenoid", "atthinmultipole", "atM66", "atQuantDiff"];

n.id = "atelemfuncs";
n.title = "Element manipulation";
n.list = ["isatelem", "atguessclass", "atshiftelem", "attiltelem"];

k.id = "atlatticefuncs";
k.title = "Lattice manipulation";
k.list=["isatlattice",...
    "-Access elements", "atindex", "atgetcells", "atgetfieldvalues", "atsetfieldvalues", ...
    "-Insert elements", "atinsertelems", "atdivelem", "atsplitelem", "insertindrift", "atsbreak",...
    "-Join elements", "atreduce", "mergedrift", "combinebypassmethod", "combinelinear45",...
    "-Other",...
    "atloadfielderrs",...
    "atsetRFCavity", "atsetshift", "atsettilt",...
    "settags", "findtags",...
    "mvelem", "mvfield", "reverse", "splitdrift"];

a.id = "atloadsave";
a.title = "Loading and Saving lattices";
a.list = ["atwritem", "atwritepy"];

b.id = "atlinearoptics";
b.title = "Linear optics";
b.list = ["atlinopt2", "atlinopt4", "atlinopt6", "beam22", "beam44", ...
    "findm44", "findm66", "findelemm44", "findelemm66"];

l.id = "atphysics";
l.title = "Physics";
l.list = ["symplectify"];

res=[m n k a b l];
end
