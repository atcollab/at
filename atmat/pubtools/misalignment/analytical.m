%% Analytical calculation of amplication factor
%
% For the quadrupole the formula follows,
%%
% 
% $$A_{x,y} = \frac{\sqrt\beta_{x,y}}{2\sqrt 2 sin(\pi \nu_{x,y})}\sqrt{\sum_i{(kl)^2_i \beta_{ix,y}}}$$
% 
aspsr_v2simple;

mach = machine_at; 
qkickx = 0; 
qkicky = 0;
% beta functions and integrated strengths at the quadrupoles
for i=1:length(THERING)
    if isfield(THERING{i},'K') & ~isfield(THERING{i},'BendingAngle')
        qkickx = qkickx + (THERING{i}.K*THERING{i}.Length).^2*mach.betax(i); 
        qkicky = qkicky + (THERING{i}.K*THERING{i}.Length).^2*mach.betay(i); 
    end
end
amp_x = sqrt(mach.betax)/(2^(3/2)*sin(pi*mach.nux(end)))*sqrt(qkickx); 
amp_y = sqrt(mach.betay)/(2^(3/2)*sin(pi*mach.nuy(end)))*sqrt(qkicky);

% Beta functions at the dipoles
betax_dip = mach.betax(findcells(THERING,'FamName','dip'));
betay_dip = mach.betay(findcells(THERING,'FamName','dip'));
dip_amp_x = THERING{15}.BendingAngle*sqrt(mach.betax)/(2^(3/2)*sin(pi*mach.nux(end)))*sqrt(sum(betax_dip)); 
dip_amp_y = THERING{15}.BendingAngle*sqrt(mach.betay)/(2^(3/2)*sin(pi*mach.nuy(end)))*sqrt(sum(betay_dip)); 

figure; 
% subplot(2,1,1,'Position',[0.15 0.25 0.7 0.7])
plot(mach.spos,abs(amp_x),'b:');
hold on;
plot(mach.spos,abs(amp_y),'r--');
a = axis;
axis([0 216/14 a(3) a(4)])
plotelementsat;
title('Quadrupole amplification factor');
xlabel('S (m)');
ylabel('Amplification \Delta x,y^{orbit}_{rms} / \Delta x,y^{magnet}_{rms}');
legend('Horizontal','Vertical');
% axes('Position',[0.1 0.2 0.7 0.7]);axis([0 216/14 a(3) a(4)])
% a = axis;
curr = get(gca,'Position');
set(gca,'Position',[curr(1) curr(2)+0.13 curr(3) curr(4)-0.13]);
a = axis;
text('Interpreter','latex',...
    'String','$$A_{x,y} = \frac{\sqrt\beta_{x,y}}{2\sqrt 2 sin(\pi \nu_{x,y})}\sqrt{\sum_i{(kl)^2_i \beta_{ix,y}}}$$',...
    'Position',[a(1)+0.25*(a(2) - a(1)), a(3)-0.25*(a(4) - a(3))]);


figure; 
% subplot(2,1,1,'Position',[0.15 0.25 0.7 0.7])
plot(mach.spos,abs(dip_amp_x),'b:');
hold on;
plot(mach.spos,abs(dip_amp_y),'r--');
a = axis;
axis([0 216/14 a(3) a(4)])
plotelementsat;
title('Dipole amplification factor');
xlabel('S (m)');
ylabel('Amplification \Delta x,y^{orbit}_{rms} / Field Error');
legend('Horizontal','Vertical');
% axes('Position',[0.1 0.2 0.7 0.7]);axis([0 216/14 a(3) a(4)])
% a = axis;
% curr = get(gca,'Position');
% set(gca,'Position',[curr(1) curr(2)+0.13 curr(3) curr(4)-0.13]);
% a = axis;
% text('Interpreter','latex',...
%     'String','$$A_{x,y} = \theta \frac{\sqrt\beta_{x,y}}{2\sqrt 2 sin(\pi \nu_{x,y})}\sqrt{\sum_i \beta_{ix,y}}}$$',...
%     'Position',[a(1)+0.25*(a(2) - a(1)), a(3)-0.25*(a(4) - a(3))]);
