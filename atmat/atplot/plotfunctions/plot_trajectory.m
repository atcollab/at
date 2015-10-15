function varargout=plot_trajectory(varargin)
%PLOT_TRAJECTORY    Plots particle trajectories
%
%Helper function for atplot: plot
%- H and V trajectories on the left axis

if nargout == 1
    nparts=size(varargin{4},2);
    refpts=1:length(varargin{2})+1;
    rout=linepass(varargin{[2 4]},refpts);
    xx=reshape(rout(1,:),nparts,[]);
    yy=reshape(rout(3,:),nparts,[]);
    labs=[arrayfun(@(n) num2str(n,'x%d'),1:nparts,'UniformOutput',false),...
        arrayfun(@(n) num2str(n,'z%d'),1:nparts,'UniformOutput',false)];
    plotdata(1).values=[xx' yy'];
    plotdata(1).labels=labs;
    plotdata(1).axislabel='x,z [m]';
    varargout={plotdata};
else
    s=findspos(varargin{1},1:length(varargin{1})+1);
    varargout={s,plot_trajectory([],varargin{:})};
end
end
