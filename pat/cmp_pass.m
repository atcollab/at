% sp3v82_lelt

%compare results of the old and new passmethods
a=ls('*pass*.c');
b=cellstr(a);
X=zeros(6,11); X(1,:) = 0.001*[-5:5];

figure(82)
for ii=1:length(b)
    indx = findcells(THERING, 'PassMethod', b{ii}(1:end-2));
    if isempty(indx)
        continue
    end
    
    Elem = THERING{indx(1)};
    fprintf('%s\n',Elem.PassMethod);
    
    cd r:\xiahuang\misc\pat\
    pX = ringpass({Elem}, X, 1); 
    
    cd R:\Controls\matlab\at\simulator\element\
    nX = ringpass({Elem}, X, 1); 

    plot(X(1,:), pX(1,:),'s-', X(1,:), nX(1,:),'o-') 
    title(Elem.PassMethod)
   
    pause
end
cd r:\xiahuang\misc\pat\


%% test tracking

>> X=zeros(6,1000); X(1,:) = 0.001*rand(1,1000);
tic; nX = ringpass(THERING, X, 1000); toc    %old
%Elapsed time is 481.753373 seconds.

>> cd r:\xiahuang\misc\pat\
>> tic; pX = ringpass(THERING, X, 1000); toc    %parallel
%Elapsed time is 116.295076 seconds.
>> norm(pX-nX)
%ans =
%     0
