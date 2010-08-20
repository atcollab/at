deletefielderr;
aspinit
numsets = 3000;

setfielderr('QFA',repmat(0,1,4),[0 1 0 0]*1e-3*1.7610967,1);
% setfielderr('QDA',repmat(0,1,4),[0 1 0 0]*1e-3*1.0715748,3);
% setfielderr('QFB',repmat(0,1,4),[0 1 0 0]*1e-3*1.5406418,3);
calcfielderr(numsets,'normal');


% ind = findcells(THERING,'FamName','QFA','QDA','QFB');
ind = findcells(THERING,'FamName','QFA');

vals = [];
for i=1:numsets
    applyfielderr quiet;
    vals(i) = getcellstruct(THERING,'PolynomB',ind(1),2);
    undofielderr quiet;
end



figure;
histvals = hist(vals,x);
hist(vals,x)
hold on;
sigma = 1e-3*1.7610967;
x = [1.755:0.0001:1.767];
px = max(histvals)*exp(-(x - 1.7610967).^2/(2*sigma^2));
plot(x,px,'o-k','LineWidth',2)
a = axis; a(1:2) = [1.75400000000000   1.77000000000000];
axis(a);
