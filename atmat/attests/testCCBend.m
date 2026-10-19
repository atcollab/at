function tests=testCCBend
%TESTCCBEND Compare with elegant ccbend1 and check the three radiation modes.
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.ref=jsondecode(fileread(fullfile(root,'pyat','test','ccbend1_reference.json')));
end

function testElegantReference(testCase)
ref=testCase.TestData.ref;p=ref.parameters;
B=atccbend(ref.family_name,p.L,p.ANGLE,p.K1,'PolynomB',[0 p.K1 p.K2/2], ...
    'fringeInt1',ref.fringeInt1','fringeInt2',ref.fringeInt2', ...
    'FringeModel',p.FRINGEMODEL,'CompensateKN',p.COMPENSATE_KN);
% Nine-point grid from ccbend.ele; inputs initially use elegant slopes.
n=-4:4;h=ref.matrix_step;q=zeros(6,1);
for j=1:6
    axis=zeros(6,9);axis(j,:)=n*h;q=[q axis]; %#ok<AGROW>
    for k=1:j-1
        pair=zeros(6,81);pair(j,:)=repelem(n,9)*h;pair(k,:)=repmat(n,1,9)*h;
        q=[q pair]; %#ok<AGROW>
    end
end
rin=q([1 2 3 4 6 5],:);norm=(1+q(6,:))./sqrt(1+q(2,:).^2+q(4,:).^2);
rin(2,:)=q(2,:).*norm;rin(4,:)=q(4,:).*norm;
out=CCBendPass(B,rin);
verifyEqual(testCase,linepass({B},rin),out,'AbsTol',1e-14);
pz=sqrt((1+out(5,:)).^2-out(2,:).^2-out(4,:).^2);
f=[out(1,:);out(2,:)./pz;out(3,:);out(4,:)./pz;out(6,:)+B.Length;out(5,:)];
w=pinv(n'.^(0:4));C=f(:,1);R=zeros(6);T=zeros(6,6,6);offset=2;
for j=1:6
    a=f(:,offset:offset+8)-C;offset=offset+9;
    R(:,j)=a*w(2,:)'/h;T(:,j,j)=a*w(3,:)'/h^2;
    for k=1:j-1
        pair=reshape(f(:,offset:offset+80)-C,6,9,9);offset=offset+81;
        for i=1:6,T(i,j,k)=w(2,:)*squeeze(pair(i,:,:))*w(2,:)'/h^2;end
    end
end
verifyEqual(testCase,C,ref.C,'AbsTol',1e-12);
verifyEqual(testCase,R,ref.R,'AbsTol',2e-11);
verifyEqual(testCase,T,ref.T,'AbsTol',2e-8);
end

function testRadiationModes(testCase)
B=atccbend('CB',1.2,.12,0,'PassMethod','CCBendRadPass','Energy',3e9);
rad=CCBendRadPass(B,zeros(6,1));verifyLessThan(testCase,rad(5),0);
diff=diffusion_matrix(B,zeros(6,1),B.Energy,atparticle('electron'),B.Length,0);
verifyGreaterThan(testCase,diff(5,5),0);
B.PassMethod='CCBendQuantPass';
quant=CCBendQuantPass(B,zeros(6,20000));
verifyEqual(testCase,var(quant(5,:),1)/diff(5,5),.973336924663,'AbsTol',.04);
verifyLessThan(testCase,abs(mean(quant(5,:))-rad(5)),5*sqrt(diff(5,5)/20000));
B.QuantumSingleParticle=0;
verifyEqual(testCase,CCBendQuantPass(B,zeros(6,1)),rad);
end
