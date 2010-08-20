%% Test Misalignment Functions
%
% The purpose is of this script is to go though a sequence of steps to
% verify that the misalignments and the statistics collected are correct.
% Begin with some initialisations.

x_error = 1e-4;
y_error = 1e-4;

%% Test CALCMISALIGN
%
% Check the distributions that are calculated by this function. Can do this
% by setting the misalignment for one element and calculating 10000 seeds
% and plot the histogram.

aspsr_v2simple;
deletemisalign;

setmisalign(7,[x_error 0 0 0 0 0],0);
calcmisalign(10000,'normal');
mis = getappdata(0,'MisalignData');
figure; hist(mis.data(7).val(1,:));
fprintf('''normal'' Mean: %f      Deviation: %f\n',...
    mean(mis.data(7).val(1,:)),std(mis.data(7).val(1,:)));
fprintf('Expect   Mean:~%f      Deviation:~%f\n',0,x_error);

calcmisalign(10000,'sho');
mis = getappdata(0,'MisalignData');
figure; hist(mis.data(7).val(1,:));
fprintf('''sho''    Mean: %f      Deviation: %f\n',...
    mean(mis.data(7).val(1,:)),std(mis.data(7).val(1,:)));

calcmisalign(10000,'equal');
mis = getappdata(0,'MisalignData');
figure; hist(mis.data(7).val(1,:));
fprintf('''equal''  Mean: %f      Deviation: %f\n',...
    mean(mis.data(7).val(1,:)),std(mis.data(7).val(1,:)));

calcmisalign(10000,'abs');
mis = getappdata(0,'MisalignData');
figure; hist(mis.data(7).val(1,:));
fprintf('''abs''    Mean: %f      Deviation: %f\n',...
    mean(mis.data(7).val(1,:)),std(mis.data(7).val(1,:)));
fprintf('Expect   Mean: %f      Deviation: %f\n',x_error,0);

%% Test APPLYMISALIGN
%
% Make sure that when the calculated misalignments are applied they are
% assigned to the correct elements and that it was done 'truthfully'. The
% first test below is just for a single element

aspsr_v2simple;
deletemisalign;

setmisalign(7,[x_error 0 0 0 0 0],0);
calcmisalign(1,'abs');
applymisalign;
fprintf('Misalign: %f   Expect: %f\n',THERING{7}.T2(1),x_error);

%% 
% Stacking misalignments
aspsr_v2simple;
deletemisalign;

setmisalign(7,[x_error 0 0 0 0 0],0);
setmisalign(7,[x_error 0 0 0 0 0],0);
setmisalign(7,[x_error 0 0 0 0 0],0);
calcmisalign(1,'abs');
applymisalign;
fprintf('Misalign: %f   Expect: %f\n',THERING{7}.T2(1),x_error*3);

%%
% Only printing out the horizontal shift. Q1 should be nonzero and all the
% rest should be 0. This test the assignment of multiple elements of one
% family.
aspsr_v2simple;
deletemisalign;

setmisalign('q1',[x_error 0 0 0 0 0],0);
calcmisalign(1,'abs');
applymisalign;
for i=1:length(THERING)
    if strcmpi(THERING{i}.FamName,'q1')
        fprintf('\n q1 (%03d): %f\n',i,THERING{i}.T2(1));
    elseif isfield(THERING{i},'T2')
        fprintf('%s(%3.2e) ',THERING{i}.FamName,THERING{i}.T2(1));
    end
end
fprintf('\nExpect horizontal misalignment of ''q1'': %f\n',x_error);
clear i;

%%
% Multiple elements of multiple families
aspsr_v2simple;
deletemisalign;

setmisalign({'q1' 'q2'},[x_error 0 0 0 0 0],0);
calcmisalign(1,'abs');
applymisalign;
for i=1:length(THERING)
    if strcmpi(THERING{i}.FamName,'q1')
        fprintf('\n q1 (%03d): %f\n',i,THERING{i}.T2(1));
    elseif strcmpi(THERING{i}.FamName,'q2')
        fprintf('\n q2 (%03d): %f\n',i,THERING{i}.T2(1));
    elseif isfield(THERING{i},'T2')
        fprintf('%s(%3.2e) ',THERING{i}.FamName,THERING{i}.T2(1));
    end
end
fprintf('\nExpect horizontal misalignment of ''q1'': %f\n',x_error);
clear i;

%%
% Prettty sure that whats output by the random number generator acutally
% gets to THERING, therefore now test something else. Calculating more than
% one set of seeds and applying each on individually and checking that
% different seeds are being loaded and that you can reload previous seeds
% that were used.

aspsr_v2simple;
deletemisalign;

setmisalign('q1',[x_error 0 0 0 0 0],0);
calcmisalign(3,'normal');
q1mis = zeros(28,3);
for j=1:3
    applymisalign;
    temp = 0;
    for i=1:length(THERING)
        if strcmpi(THERING{i}.FamName,'q1')
            temp = temp + 1;
            q1mis(temp,j) = THERING{i}.T2(1);
        end
    end
end

% Calculate the correlation factor for the three different combinations to
% determine if any of the sets of seeds are repeated. Correlation factor
% should be small for all three
disp('Expect samll correlation factors, k, for the three cases:');
k = mean(q1mis(:,1).*q1mis(:,2))/...
    sqrt(mean(q1mis(:,1).^2)*mean(q1mis(:,2).^2))
k = mean(q1mis(:,2).*q1mis(:,3))/...
    sqrt(mean(q1mis(:,2).^2)*mean(q1mis(:,3).^2))
k = mean(q1mis(:,1).*q1mis(:,3))/...
    sqrt(mean(q1mis(:,1).^2)*mean(q1mis(:,3).^2))

% Apply the seed number 2 again and store the misalignments from THERING.
% Comparison of the current misalignment and that stored in q1mis(:,2)
% should produce a correlation factor of 1.
applymisalign(2);
temp = 0;
for i=1:length(THERING)
    if strcmpi(THERING{i}.FamName,'q1')
        temp = temp + 1;
        q1mis(temp,4) = THERING{i}.T2(1);
    end
end
disp('Expect correlation of 1 below:')
k = mean(q1mis(:,2).*q1mis(:,4))/...
    sqrt(mean(q1mis(:,2).^2)*mean(q1mis(:,4).^2))

clear q1mis k temp j;

%% Test PRINTMISALIGN
%
% The function prints out to the screen or a file the misalignments of each
% individual misalignment or collects the mean and deviation statistics for
% a family of elements. This can be verified by setting the misalignments
% all to one value only, therefore we know the mean and deviation (which
% will be 0). Secondly use a normal distribution with a sigma of X, and the
% statistics should also show a deviation ~X.
aspsr_v2simple;
deletemisalign;

setmisalign('q1',[x_error 0 0 0 0 0],0);
calcmisalign(1,'abs');
applymisalign;
printmisalign q1
disp(['Expect to see an x deviation of zero and a mean of ' num2str(x_error)])

%%
% Display the individual magnet misalignments
printmisalign q1 ind

%%
% Display the individual magnet misalignments for all elements
printmisalign ind file test.txt
edit test.txt

%%
% Calculate the statistics for the normal distribution.
calcmisalign(1,'normal');
applymisalign;
printmisalign q1
disp(['Expect to see an x deviation around zero and a mean of ~' num2str(x_error)])


%% Test Misalignment Operations
%
% By now we should be confident that the various tools put in place
% faithfully does the translations, however there are rotations that have
% to be checked out to ensure that they follow the right hand rule for
% rotation, especially for the girders. 
%
% Additionally some tests will be carried out to ensure that the shifts and
% rotations are realistic by looking at the effect on the optics.
%
%% Test SETMISALIGN rotations in x and y
%
% A rotation about x axis will result in an effective translation of the
% beam coordinates and 'angle' in y, to first order. Assuming that these
% rotations are small then the change in 's' will be small. Which is why
% the R1 and R2 matrices are not being used for this purpose.

aspsr_v2simple;
deletemisalign;

xshift = 1e-3;
yshift = 2e-3;
setmisalign(7,[xshift 0 yshift 0 0 0],0);
calcmisalign(1,'abs')
applymisalign
fprintf('T1 = [%5.3e %5.3e %5.3e %5.3e %5.3e %5.3e]\n',THERING{7}.T1);
fprintf('T2 = [%5.3e %5.3e %5.3e %5.3e %5.3e %5.3e]\n',THERING{7}.T2);
fprintf('Expect:\n');
fprintf('T1 = [%5.3e 0 %5.3e 0 0 0]\n',-xshift,-yshift);
fprintf('T2 = [%5.3e 0 %5.3e 0 0 0]\n',xshift,yshift);


%%
aspsr_v2simple;
deletemisalign;

xrot = 3e-3;
yrot = 3e-3;
setmisalign(7,[0 xrot 0 0 0 0],0);
calcmisalign(1,'abs')
applymisalign
fprintf('T1 = [%5.3e %5.3e %5.3e %5.3e %5.3e %5.3e]\n',THERING{7}.T1);
fprintf('T2 = [%5.3e %5.3e %5.3e %5.3e %5.3e %5.3e]\n',THERING{7}.T2);
fprintf('Expect:\n');
fprintf('T1 = [0 0 %5.3e %5.3e 0 0]\n',-THERING{7}.Length*0.5*tan(xrot),xrot);
fprintf('T2 = [0 0 %5.3e %5.3e 0 0]\n',THERING{7}.Length*0.5*tan(xrot),-xrot);

deletemisalign;
setmisalign(7,[0 0 0 yrot 0 0],0);
calcmisalign(1,'abs')
applymisalign
fprintf('T1 = [%5.3e %5.3e %5.3e %5.3e %5.3e %5.3e]\n',THERING{7}.T1);
fprintf('T2 = [%5.3e %5.3e %5.3e %5.3e %5.3e %5.3e]\n',THERING{7}.T2);
fprintf('Expect:\n');
fprintf('T1 = [%5.3e %5.3e 0 0 0 0]\n',THERING{7}.Length*0.5*tan(yrot),-yrot);
fprintf('T2 = [%5.3e %5.3e 0 0 0 0]\n',-THERING{7}.Length*0.5*tan(yrot),yrot);

%%
aspsr_v2simple;
srot = pi/rand;
deletemisalign;
setmisalign(7,[0 0 0 0 0 srot],0);
calcmisalign(1,'abs')
applymisalign
fprintf('%5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n',THERING{7}.R1);
fprintf('Expect\n');
C = cos(-srot); S = sin(-srot);
temp = [C  0 -S  0  0  0;
        0  C  0 -S  0  0;
        S  0  C  0  0  0;
        0  S  0  C  0  0;
        0  0  0  0  1  0;
        0  0  0  0  0  1];
fprintf('%5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n',temp);


%% Test SETGIRDERMISALIGN
%
% Start with a simple transverse shift, then examine the girder rotations
% about x and y axes. The s roataion is trivial.

aspsr_v2simple;
deletemisalign;

xshift = 1e-7;
yshift = 2e-7;
setgirdermisalign('g1m1','g1m2',[xshift 0 yshift 0 0 0],0);
calcmisalign(1,'normal')
applymisalign
printmisalign q1 s1 s2 ind
printmisalign q1 s1 s2

%%
% X rotation of the girder.
aspsr_v2simple;
deletemisalign;

xrot = 1e-3;
setgirdermisalign('g1m1','g1m2',[0 xrot 0 0 0 0],0);
calcmisalign(1,'abs')
applymisalign
printmisalign q1 s1 s2 ind
s11pos = 2.6983+0.05;
s12pos = 2.7983+0.05;
q1pos = 3.0883+.355/2;
s21pos = 3.6083+0.05;
s22pos = 3.7083+0.05;
girdercentre = (3.7083+0.1 + 2.6983)/2;
fprintf('Expect: %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e \n',...
    0,xrot,(girdercentre - s11pos)*sin(xrot),0,(girdercentre - s11pos)*(1-cos(xrot)),0);
fprintf('Expect: %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e \n',...
    0,xrot,(girdercentre - s12pos)*sin(xrot),0,(girdercentre - s12pos)*(1-cos(xrot)),0);
fprintf('Expect: %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e \n',...
    0,xrot,(girdercentre -  q1pos)*sin(xrot),0,(girdercentre -  q1pos)*(1-cos(xrot)),0);
fprintf('Expect: %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e \n',...
    0,xrot,(girdercentre - s21pos)*sin(xrot),0,(girdercentre - s21pos)*(1-cos(xrot)),0);
fprintf('Expect: %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e \n',...
    0,xrot,(girdercentre - s22pos)*sin(xrot),0,(girdercentre - s22pos)*(1-cos(xrot)),0);

%%
% Y rotation of the girder.
aspsr_v2simple;
deletemisalign;

yrot = 2e-3;
setgirdermisalign('g1m1','g1m2',[0 0 0 yrot 0 0],0);
calcmisalign(1,'abs')
applymisalign
printmisalign q1 s1 s2 ind
s11pos = 2.6983+0.05;
s12pos = 2.7983+0.05;
q1pos = 3.0883+.355/2;
s21pos = 3.6083+0.05;
s22pos = 3.7083+0.05;
girdercentre = (3.7083+0.1 + 2.6983)/2;
fprintf('Expect:\n');
fprintf('s1: %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e \n',...
    (girdercentre - s11pos)*sin(-yrot),0,0,yrot,(girdercentre - s11pos)*(1-cos(-yrot)),0);
fprintf('s1: %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e \n',...
    (girdercentre - s12pos)*sin(-yrot),0,0,yrot,(girdercentre - s12pos)*(1-cos(-yrot)),0);
fprintf('q1: %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e \n',...
    (girdercentre -  q1pos)*sin(-yrot),0,0,yrot,(girdercentre -  q1pos)*(1-cos(-yrot)),0);
fprintf('s2: %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e \n',...
    (girdercentre - s21pos)*sin(-yrot),0,0,yrot,(girdercentre - s21pos)*(1-cos(-yrot)),0);
fprintf('s2: %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e  |  %5.3e \n',...
    (girdercentre - s22pos)*sin(-yrot),0,0,yrot,(girdercentre - s22pos)*(1-cos(-yrot)),0);

%%
% S rotation of the girder.
