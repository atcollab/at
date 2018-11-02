function varargout = pyproxy(funcname,ring,varargin)
%Convert structures to a form accessible in python

func=str2func(funcname);
[outp{1:nargout}]=func(ring(:), varargin{:});
varargout=cellfun(@convert,outp,'UniformOutput',false);

    function argout=convert(argin)
        if isstruct(argin)
            argout=struct();
            fnames=fieldnames(argin);
            for i=1:length(fnames)
                fn=fnames{i};
                % pack values along the 3rd dimension and move refpoints to
                % 1st dimension
                argout(1).(fn)=permute(cat(3,argin.(fn)),[3 1 2]);
            end
        else
            argout=argin;
        end
    end
end

