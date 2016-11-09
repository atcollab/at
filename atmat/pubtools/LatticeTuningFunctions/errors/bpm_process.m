function bpmreading = bpm_process(bpmorbit,rel,tel,trand)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% by L.F. Jun 2016, ESRF,  K:\machine\matlab\atlf

bp=cellfun(@process,mat2cell(bpmorbit,2,ones(1,size(bpmorbit,2))),rel,tel,trand,'UniformOutput',false);
bpmreading=cat(2,bp{:});

    function x=process(x0,rel,tel,trand)
        x=rel*x0 + tel + trand.*randn(2,1);
    end
end
