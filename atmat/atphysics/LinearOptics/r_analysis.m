function [bk0,phi,bks,as]=r_analysis(a0,ms)

nv=size(a0,1);
dms=nv/2;
slices=num2cell(reshape(1:nv,2,dms),1);
s=jmat(dms);

astd=standardize(a0);

bk0=cellfun(@(slc) astd(:,slc)*astd(:,slc)',slices,'UniformOutput',false);
bk0=cat(3,bk0{:});
[phi,bks,as]=cellfun(@(mi) propagate(mi*astd),ms,'UniformOutput',false);

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
            rot=-getphase(aa(range,range));
            cs=cos(rot);
            sn=sin(rot);
            R(range,range)=[cs sn;-sn cs];
        end
    end

    function [phi,bk,ai]=propagate(ai)
        % Propagate the A matrices
        ais=ai*s;
        invai=ai\s';
        bk=cellfun(@(slc) ais(:,slc)*invai(slc,:),slices,'UniformOutput',false);
%         bk=cellfun(@(slc) ai(:,slc)*ai(:,slc)',slices,'UniformOutput',false);
        bk=cat(3,bk{:});
        phi=cellfun(@(slc) getphase(ai(slc,slc)),slices);
    end
end

