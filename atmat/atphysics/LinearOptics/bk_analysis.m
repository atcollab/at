function varargout=bk_analysis(m66,ms)

nv=size(m66,1);
dms=nv/2;
slices=num2cell(reshape(1:nv,2,dms),1);
Tk=cellfun(@makeTk,slices,'UniformOutput',false);

astd=standardize(amat(m66));

mms=squeeze(num2cell(ms,[1 2]));

[varargout{1:dms+1}]=cellfun(@(mi) propagate(mi*astd),mms,'UniformOutput',false);

    function tk=makeTk(slc)
        % Build the Tk selection matrices
        tk=zeros(nv,nv);
        tk(slc,slc)=eye(2);
    end

    function phi=getphase(a22)
        % Return the phase for A standardization
        phi=atan2(a22(1,2),a22(1,1));
    end

    function RR=standardize(aa)
        % Apply rotation to set A in std form
        R=zeros(nv,nv);
        cellfun(@rot22,slices);
        RR=aa*R;

        function rot22(range)
            r1=-getphase(aa(range,range));
            cs=cos(r1);
            sn=sin(r1);
            R(range,range)=[cs sn;-sn cs];
        end
    end

    function varargout=propagate(ai)
        % Propagate the A matrices
        bk=cellfun(@(ik) ai*ik*ai',Tk,'UniformOutput',false);
        phi=cellfun(@(slc) getphase(ai(slc,slc)),slices);
        varargout=[{phi} bk];
    end
end

