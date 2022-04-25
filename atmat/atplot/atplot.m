function varargout = atplot(varargin)
%ATPLOT Plots optical functions
%
%ATPLOT                 Plots THERING in the current axes
%
%ATPLOT(RING)           Plots the lattice specified by RING
%
%ATPLOT(AX,RING)        Plots in the axes specified by AX
%
%ATPLOT(AX,RING,DPP)    Plots at momentum deviation DPP
%
%ATPLOT(...,[SMIN SMAX])  Zoom on the specified range
%
%ATPLOT(...,'OptionName',OptionValue,...) Available options:
%   'twiss_in',structure of optics   transferline optics
%               (ex: [optics_struct,~,~]=atlinopt(ring,0,1);
%              atplot(ring,[0 10],@plBeamSize,'twiss_in',optics_struct);)
%   'comment',true|false        Prints lattice information (default:true)
%   'synopt',true|false         Plots the lattice elements
%   'labels',REFPTS             Display the names of selected element names
%   'index',REFPTS             Display the index of selected element names
%   'leftargs',{properties}     properties set on the left axis
%   'rightargs',{properties}    properties set on the right axis
%   'KeepAxis'                  flag to keep R1,R2,T1,T2 at each slice in
%                               detailed plots (mandatory with vert. bend).
%
%ATPLOT(...,@PLOTFUNCTION,args...)
%	Allows for a user supplied function providing the values to be plotted
%   PLOTFUNCTION must be of form:
%   PLOTDATA=PLOTFUNCTION(LINDATA,RING,DPP,args...)
%
%       PLOTDATA: structure array,
%         PLOTDATA(1) describes data for the left (main) axis
%           PLOTDATA(1).VALUES: data to be plotted, length(ring)+1 x nbcurve
%           PLOTDATA(1).LABELS: curve labels, cell array, 1 x nbcurve
%           PLOTDATA(1).AXISLABEL: string
%         PLOTDATA(2) optional, describes data for the right (secondary) axis
%
%   The default function displayed below as an example plots beta functions
%	and dispersion
%
% function plotdata=defaultplot(lindata,ring,dpp,varargin)
% beta=cat(1,lindata.beta);                     % left axis
% plotdata(1).values=beta;
% plotdata(1).labels={'\beta_x','\beta_z'};
% plotdata(1).axislabel='\beta [m]';
% dispersion=cat(2,lindata.Dispersion)';        % right axis
% plotdata(2).values=dispersion(:,1);
% plotdata(2).labels={'\eta_x'};
% plotdata(2).axislabel='dispersion [m]';
% end
%
%CURVE=ATPLOT(...) Returns handles to some objects:
%   CURVE.PERIODICIY ring periodicity
%   CURVE.LENGTH    structure length
%   CURVE.DPP       deltap/p
%   CURVE.LEFT      Handles to the left axis plots
%   CURVE.RIGHT     Handles to the right axis plots
%   CURVE.LATTICE   Handles to the Element patches: structure with fields
%          Dipole,Quadrupole,Sextupole,Multipole,BPM,Label
%   CURVE.COMMENT   Handles to the Comment text
%
%
%ATPLOT calls the more general ATBASEPLOT function, which uses a slightly
%different syntax.
%
%See also: atbaseplot

global THERING %#ok<GVMIS>

% Select axes for the plot
[ax,args]=getargs(varargin,[],'check',@(x) isscalar(x) && ishandle(x));
if isempty(ax), ax=gca; end

% Select the lattice
[ring,args]=getargs(args,THERING,'check',@iscell);

% Extract plotfunction and plotargs
funcarg=find(cellfun(@(arg) isa(arg,'function_handle'),args),1);
if isempty(funcarg)
    funcarg=length(args)+1;
    options={@plotbetadisp};
else
    options=args(funcarg:end);
end
args=getdparg(args(1:funcarg-1));
[varargout{1:nargout}] = wrapper6d(ring,@doplot,ax,options, args{:});


    function varargout=doplot(ring,is6d,ax,options,varargin)
        [comment,varargs]=getoption(varargin,'comment',true);
        varargs=cellfun(@correct, varargs, 'UniformOutput',false);
        [dpargs,varargs]=getoption(varargs,{'dp','dct','df'});
        [srange,varargs]=getargs(varargs,[0 inf],'check',@(x) isnumeric(x) && numel(x)==2);
        [largs,varargs]=opticsoptions(varargs);

        lindata=[];

        curve=xplot(ring,is6d,ax,srange,@ringplot,[options largs],dpargs{:},varargs{:});

        if comment && ~isempty(curve.left)
            set(ax,'Position',[.13 .11 .75 .75]);
            tuneper=lindata(end).mu/2/pi;
            tunes=curve.periodicity*tuneper;
            circumference=curve.periodicity*curve.length;
            if curve.periodicity > 1, plural='s'; else plural=''; end %#ok<SEPEX>
            if isfinite(curve.dpp)
                dppstring=sprintf('\\deltap/p=%.3f',curve.dpp);
            else
                dppstring='';
            end
            line1=sprintf('\\nu_x=%8.3f      %s',tunes(1),dppstring);
            line2=sprintf('\\nu_z=%8.3f      %2i %s, C=%10.3f',tunes(2),...
                curve.periodicity,['period' plural],circumference);
            curve.comment=text(ax,-0.14,1.12,{line1;line2},'Units','normalized',...
                'VerticalAlignment','top');
        end
        if nargout>0, varargout={curve}; end

        function [s,plotdata]=ringplot(ring,dpp,plotfun,varargin)
            [linargs,vargs] = opticsoptions(varargin);
            [ringdata,lindata]=atlinopt6(ring,1:length(ring)+1,linargs{:}); %#ok<ASGLU>
            s=cat(1,lindata.SPos);
            plotdata=plotfun(lindata,ring,dpp,vargs{:});
        end

        function opt=correct(opt)
            % Replace the obsolete 'inputtwiss' keyword with 'twiss_in'
            if strcmp(opt, 'inputtwiss')
                opt='twiss_in';
            end
        end
    end

end
