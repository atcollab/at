testline dipole;

figure;
td = linepass(THERING,[0 0 0 0 0 0]',[1:length(THERING)]);
plot(td(3,:))
setsrot(2,pi/4);
td = linepass(THERING,[0 0 0 0 0 0]',[1:length(THERING)]);
hold on;
plot(td(3,:))



