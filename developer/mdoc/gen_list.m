function [dd,pp] = gen_list(d)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

p=[];
[~,dd,~]=fileparts(d);
files = dir(fullfile(d,'*.m'))';
fnames = arrayfun(@(a) string(fullfile(a.folder,a.name)),files);
for name=fnames
    try
        hnext=h1_line(name);
        p=[p hnext]; %#ok<AGROW>
    catch err
        disp(err.message)
    end
end
pp=[p.name];
disp("["""+strjoin(pp,'", "')+"""]");
end
