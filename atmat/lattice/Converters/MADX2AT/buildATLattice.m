function Lat=buildATLattice(list,s,totL)
% given a list (cell array) of elements with specified field Spos (center of element (madx default)) in a
% vector returns a cell array with elements without Spos field and
% appropriate Drifts spaces between. Drifts of the same length have the same name.
Lat={};
latind=1;
zero=0;

[s,indpos]=sort(s);
list=list(indpos);
elemalign='b'; % (c) center (b) begin (e) end

for i=1:length(list)
    L=list{i}.('Length');
    
 
    spos=s(i)-L/2; % refer to element entrance
    D=[];% drift name
    
    switch elemalign
        case 'c'
            DL=spos-zero-L/2; % drift length
        case 'b'
            DL=spos-zero; % drift length
        case 'e'
            DL=spos-zero-L; % drift length
        otherwise
            DL=spos-zero-L/2; % drift length
    end

    
    if DL>1e-7
        Dnam=['DR_' num2str(DL)]; %['Drift'];
        Dnam(ismember(Dnam,'.'))=[];
        D=atdrift(Dnam,DL);
        Lat{latind}=D;
        Lat{latind+1}=list{i};
        latind=latind+1;
        
    elseif (DL>=0 && DL<1e-7)
        Lat{latind}=list{i};
        
    elseif DL>-1e-5 % fix roundings
      
        Lat{latind}=list{i};
        Lat{latind-1}.('Length')=Lat{latind-1}.('Length')+DL;
    
    else
      % list{i}
        disp([ 'negative drift: ' num2str(DL)])
    end
    



    latind=latind+1;
    
    switch elemalign
        case 'c'
            zero=spos+L/2;
        case 'b'
            zero=spos+L;
        case 'e'
            zero=spos;
        otherwise
            zero=spos+L/2;
    end
    
   
    
    
end


% add last dirft to reach lattice length

Dnam=['DR_' num2str(totL-zero)];%['Drift'];
        Dnam(ismember(Dnam,'.'))=[];
D=atdrift(Dnam,totL-zero);
       
Lat{latind}=D;


Lat=Lat'; % return column

% clean unusefull drifts

%
%LD=getcellstruct(Lat,'Length',drifts);

%Lat(drifts(LD<1e-7))=[];






