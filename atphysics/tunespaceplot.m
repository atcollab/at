function R = tunespaceplot(XTUNE,YTUNE,RESORDER,varargin);
%TUNESPACEPLOT draws a tune diagram
% resonance lines: m*nu_x + n*nu_y = p

%
% TUNESPACEPLOT(XTUNE, YTUNE, ORDER)
%
% XTUNE = [XTUNEMIN,XTUNEMAX] 
% YTUNE = [YTUNEMIN,YTUNEMAX] - plotting range in tune space
% RESORDER - resonance order: |m| + |n|

% TUNESPACEPLOT(XTUNE, YTUNE, ORDER, FIGHANDLE)
% TUNESPACEPLOT(XTUNE, YTUNE, ORDER, FIGHANDLE)


% 


 

if nargin>3
    if ishandle(varargin{1}) & strcmp(get(varargin{1},'Type'),'figure')
        % Plot tune space plot
        figure(varargin{1});
    else % create new figure
        figure
        axes;   
    end
end
if nargin>4
    LINEARGS = varargin(2:end);
else
    LINEARGS = {};
end


axis([XTUNE,YTUNE]);
axis square
        

R = zeros(8*length(RESORDER),6);
NLMAX = 0;
for r = RESORDER
    for m = 0:r
        n = r-m;
        
        % Lower
        p1 = ceil(m*XTUNE(1)+n*YTUNE(1));
        p2 = floor(m*XTUNE(2)+n*YTUNE(1));
        
            
        for p =p1:p2 
            if m % lines with m=0 do not cross upper and lower sides 
                NLMAX = NLMAX+1;
                R(NLMAX,:) = [abs(m)+abs(n),m,n,p,(p-n*YTUNE(1))/m,YTUNE(1)];
            end
        end
        
        % Left
        p1 = ceil(m*XTUNE(1)+n*YTUNE(1));
        p2 = floor(m*XTUNE(1)+n*YTUNE(2));
        
        
        for p =p1:p2
            if n % lines with n=0 do not cross left and right sides 
                NLMAX = NLMAX+1;
                R(NLMAX,:) = [abs(m)+abs(n),m,n,p,XTUNE(1),(p-m*XTUNE(1))/n];
            end
        end
        
        % Upper
        p1 = ceil(m*XTUNE(1)+n*YTUNE(2));
        p2 = floor(m*XTUNE(2)+n*YTUNE(2));
        
        for p=p1:p2
            if m
                NLMAX = NLMAX+1;
                R(NLMAX,:) = [abs(m)+abs(n),m,n,p,(p-n*YTUNE(2))/m,YTUNE(2)];
            end
        end
        
        % Right
        p1 = ceil(m*XTUNE(2)+n*YTUNE(1));
        p2 = floor(m*XTUNE(2)+n*YTUNE(2));
        
        for p=p1:p2
            if n
                NLMAX = NLMAX+1;
                R(NLMAX,:) = [abs(m)+abs(n),m,n,p,XTUNE(2),(p-m*XTUNE(2))/n];
            end
        end
        
        % ======================== 
        n = -r+m;
        
        % Lower
        p1 = ceil(m*XTUNE(1)+n*YTUNE(1));
        p2 = floor(m*XTUNE(2)+n*YTUNE(1));
        
        for p =p1:p2 
            if m % lines with m=0 do not cross upper and lower sides 
                NLMAX = NLMAX+1;
                R(NLMAX,:) = [abs(m)+abs(n),m,n,p,(p-n*YTUNE(1))/m,YTUNE(1)];
            end
        end
        
        % Left
        % Note: negative n
        p1 = floor(m*XTUNE(1)+n*YTUNE(1));
        p2 = ceil(m*XTUNE(1)+n*YTUNE(2));
        
        for p =p2:p1
            if n % lines with n=0 do not cross left and right sides 
                NLMAX = NLMAX+1;
                R(NLMAX,:) = [abs(m)+abs(n),m,n,p,XTUNE(1),(p-m*XTUNE(1))/n];
            end
        end
        
        % Upper
        p1 = ceil(m*XTUNE(1)+n*YTUNE(2));
        p2 = floor(m*XTUNE(2)+n*YTUNE(2));
        
        for p=p1:p2
            if m
                NLMAX = NLMAX+1;
                R(NLMAX,:) = [abs(m)+abs(n),m,n,p,(p-n*YTUNE(2))/m,YTUNE(2)];
            end
        end
        
        % Right
        % Note: negative n
        
        p1 = floor(m*XTUNE(2)+n*YTUNE(1));
        p2 = ceil(m*XTUNE(2)+n*YTUNE(2));
        for p=p2:p1
            if n
                NLMAX = NLMAX+1;
                R(NLMAX,:) = [abs(m)+abs(n),m,n,p,XTUNE(2),(p-m*XTUNE(2))/n];
            end
        end
    end
end
%R = sortrows(R(1:NLMAX,:));
R = unique(R(1:NLMAX,:),'rows');
[temp,I,J] = unique(R(:,1:4),'rows');
K = I(find(diff([0;I])==2))-1;

RESNUM = [R(K,1:4)]; % [o, m, n, p] O = |m| + |n|
X1 = R(K,5);
X2 = R(K+1,5);
Y1 = R(K,6);
Y2 = R(K+1,6);


% Remove accidental lines that are on the box edge
K1 = (X1==X2) & (X1==XTUNE(1));
K2 = (X1==X2) & (X1==XTUNE(2));
K3 = (Y1==Y2) & (Y1==YTUNE(1));
K4 = (Y1==Y2) & (Y1==YTUNE(2));

K = find(~(K1 | K2 | K3 | K4));


RESNUM = RESNUM(K,:);
X1 = X1(K);
X2 = X2(K);
Y1 = Y1(K);
Y2 = Y2(K);


R = [RESNUM,X1,Y1,X2,Y2];






NL = size(RESNUM,1);
for i = 1:NL
    hl = line([X1(i) X2(i)],[Y1(i) Y2(i)]);
    if ~isempty(LINEARGS)
        set(hl,LINEARGS{:});
    end
end


    