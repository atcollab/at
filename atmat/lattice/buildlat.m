function buildlat(ELIST)
%BUILDLAT places elements from FAMLIST into cell array THERING
% in the order given by ineger arry ELIST
% to be use in Accelerator Toolbox lattice definition files

global FAMLIST THERING
THERING=cell(size(ELIST(:)));
for i=1:length(THERING)
   THERING{i} = FAMLIST{ELIST(i)}.ElemData;
   FAMLIST{ELIST(i)}.NumKids=FAMLIST{ELIST(i)}.NumKids+1;
   FAMLIST{ELIST(i)}.KidsList = [FAMLIST{ELIST(i)}.KidsList i];
end
