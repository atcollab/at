function [phi,bks,as]=r_analysis(a0,ms)

nv=size(a0,1);
dms=nv/2;
slices=num2cell(reshape(1:nv,2,dms),1);
s=kmat(dms);

astd=standardize(a0);

sw=warning;
warning('off','MATLAB:illConditionedMatrix');
[phi,bks,as]=cellfun(@(mi) propagate(mi*astd),ms,'UniformOutput',false);
warning(sw);

    function phi=getphase(a22)
        % Return the phase for A standardization
        phi=atan2(a22(1,2),a22(1,1));
    end

    function astd=standardize(a0)
        % Apply rotation to set A in std form
        allv=cellfun(@rot22,slices,'UniformOutput',false);
        astd=cat(2,allv{:});

        function v=rot22(range)
            rot=-getphase(a0(range,range));
            cs=cos(rot);
            sn=sin(rot);
            v=a0(:,range)*[cs sn;-sn cs];
        end
    end

    function [phi,bk,ai]=propagate(ai)
        % Propagate the A matrices
        ais=ai*s;
        invai=ai\s';
        bk=cellfun(@(slc) ais(:,slc)*invai(slc,:),slices,'UniformOutput',false);
%       bk=cellfun(@(slc) ai(:,slc)*ai(:,slc)',slices,'UniformOutput',false);   % Only if symplectic
        bk=cat(3,bk{:});
        phi=cellfun(@(slc) getphase(ai(slc,slc)),slices);
    end

    function mat=kmat(dim)
        S2 = [0 1; -1 0];
        
        if(dim==1)
            mat=S2;
        elseif(dim==2)
            mat = blkdiag(S2,S2);
        elseif(dim==3)
            mat = blkdiag(S2,S2,S2');
        else
            error('Dim is 1, 2 or 3')
        end
    end
end
