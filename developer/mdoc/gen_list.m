function pp = gen_list(pth)
%GEN_LIST	Display a list of still non-documented AT functions

chapters=atchapters();
known=[chapters.contents];
p=[];
folders=split(string(genpath(pth)),pathsep)';
for d=folders(1:end-1)
    files = dir(fullfile(d,'*.m'))';
    fnames = arrayfun(@(a) string(fullfile(a.folder,a.name)),files);
    for name=fnames
        [~,nm,~]=fileparts(name);
        if ~any(strcmp(nm,known))
            try
                hnext=h1_line(name);
                p=[p hnext]; %#ok<AGROW>
            catch err
                disp(err.message)
            end
        end
    end
end
pp=[p.name];
disp("["""+strjoin(pp,'", "')+"""]");
end
