function res = atchapters()
%ATCHAPTERS     Describe the User Guide chapters
%
%CHAPTERS=ATCHAPTERS()
%
%CHAPTERS is a structure array where each element defines one chapter
%CHAPTERS has 3 fields;
%   id:         identifiers used as file name for HTML files
%   title:      Chapter title
%   contents:   String array. Each string defines an item with the
%               following convention:
%                   string starting with "-": section header
%                   string starting with "0": plain text
%                   otherwise: name of an AT function

m.id = "atelemcreate";
m.title = "Element creation";
m.contents = ["atbaselem", "atelem", "atringparam", "atdrift", "atquadrupole",...
    "atrbend", "atsbend", "atsextupole", "atskewquad", "atmultipole", "atrfcavity",...
    "atmarker", "atmonitor", "ataperture", "atcorrector", "atidtable", "atwiggler",...
    "atdampMatElem", "atsolenoid", "atthinmultipole", "atM66", "atQuantDiff"];

n.id = "atelemfuncs";
n.title = "Element manipulation";
n.contents = ["isatelem", "atguessclass", "atshiftelem", "attiltelem"];

k.id = "atlatticefuncs";
k.title = "Lattice manipulation";
k.contents=["isatlattice",...
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
a.contents = ["-Binary files", ...
    "0Lattices can be saved as binary mat-files using the standard load and save commands", ...
    "-Text Files", "atwritem", "atwritepy"];

b.id = "atlinearoptics";
b.title = "Linear optics";
b.contents = ["-Closed orbit" "findorbit", "findorbit4", "findorbit6", "findsyncorbit", ...
    "-Transfer matrices", "findm44", "findm66", "findelemm44", "findelemm66", ...
    "-Optical functions", "atlinopt2", "atlinopt4", "atlinopt6", "beam22", "beam44"];

c.id = "atradiation";
c.title = "Radiation";
c.contents = ["check_radiation", "atenergy", "atgetU0", "atdampingrates", "atradon", "atradoff", ...
    "quantumDiff", "ohmienvelope", "DipoleRadiation", "WigglerRadiation"];

d.id = "atsummary";
d.title = "Parameter summary";
d.contents = ["atx", "atsummary", "ringpara"];


l.id = "atphysics";
l.title = "Physics";
l.contents = ["symplectify"];

res=[m n k a b c d l];
end
