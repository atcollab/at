function C = findrespmatmisal(RING, OBSINDEX, PERTURB, PVALUE, field1,m1,n1,misal,inCOD)
%
% 
% See also ATPARAMGROUP, FINDORBIT, FINDORBIT4, FINDORBIT6, FINDSYNCORBIT

warning('off','all');

O = length(OBSINDEX);
P = length(PERTURB);
C = {zeros(O,P),zeros(O,P),zeros(O,P),zeros(O,P)};

if length(PVALUE) ~= P
    PVALUE = PVALUE(ones(1,P(1)));
end


    if nargin < 9
        error('Incorrect number of inputs');
    end
    
    if ~ischar(field1) % Check that the FIELD argument is a string
        error('The 5-th argument FIELD must be a string');
    end
    
    if ~isnumeric(m1) | length(m1)>1 % Check that the M argument is a scalar
        error('The 6-th argument FIELD must be a scalar');
    end
    if ~isnumeric(n1) | length(n1)>1 % Check that the M argument is a scalar
        error('The 6-th argument FIELD must be a scalar');
    end


    orbit_function_handle = @findorbit6Err;
    orbit_function_args   = {OBSINDEX,inCOD};
    
    
    %ORBIT = findorbit4(RING,0,OBSINDEX);
    
    ORBIT = feval(orbit_function_handle,RING,orbit_function_args{:});
    
    mn = {m1,n1};
    dx=misal;
    dy=misal;
    zz=0;
    
    for i = 1:P
        
        oldvalue = atgetfieldvalues(RING,PERTURB(i),field1,mn);
        RINGx = atsetfieldvalues(RING,PERTURB(i),field1,mn,oldvalue+PVALUE(i));
        RINGx = atsetshift(RINGx,PERTURB(i),dx,zz);
        ORBITx  = feval(orbit_function_handle,RINGx,orbit_function_args{:});
        DORBITx = (ORBITx - ORBIT);
       
        oldvalue = atgetfieldvalues(RING,PERTURB(i),field1,mn);
        RINGy = atsetfieldvalues(RING,PERTURB(i),field1,mn,oldvalue+PVALUE(i));
        RINGy = atsetshift(RINGy,PERTURB(i),zz,dy);
        ORBITy  = feval(orbit_function_handle,RINGy,orbit_function_args{:});
        DORBITy = (ORBITy - ORBIT);
       
        C{1}(:,i) = DORBITx(1,:);
        C{2}(:,i) = DORBITx(2,:);
        C{3}(:,i) = DORBITy(3,:);
        C{4}(:,i) = DORBITy(4,:);
    end

warning('on','all');

