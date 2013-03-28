function atDisplayConstraintsChange(Constraints,results)
% This funciton evaluates the contraints defined in Constraints for lattice
% Constraints: cell array of struct('Fun',@functname,'Min',min,'Max',max,'OtherParameters',otherargs}
%
% Constraints: structure array struct('Fun',@(ring)functname(ring,parameters),
%                                   'Min',min, % of unweigthed funct val
%                                   'Max',max,
%                                   'Weight',w,
%                                   'RefPoints',[])
%
% Constraints: structure array struct(...
%                     'Fun',@(ring,lindata,globaldata,refpts)functname(...
%                             ring,lindata,globaldata,refpts,parameters),
%                                   'Min',min, % of unweigthed funct val
%                                   'Max',max,
%                                   'Weight',w,
%                                   'RefPoints',refpts);
%
% lindata is the output of atlinopt
% globdata.tune=tune fromk atlinopt
% globdata.chrom=chrom from atlinopt
%
% functname: handle to vector valued function: [res]=functname(THERING,otherargs)
%
% min and max have to be the same size as res. (this allows to give functions as limits!)
%
% created 30-8-2012
% updated 12-10-2012 other_function_args is a cell array, ifit is not it is
%                    transformed in a cell array


disp('   ')
disp('Final constraints values:')
disp('   ')
disp('Name          lat_indx      before         after           low            high       min dif before    min dif after  ')
ok=arrayfun(@dispc, Constraints,results); %#ok<NASGU>
disp('   ')
disp('-----oooooo----oooooo----oooooo----')
disp('    ')

    function ok=dispc(cstr,res)
        nv=length(res.val1);
        ConstrName=sprintf('%-20.20s',func2str(cstr.Fun));
        disp([repmat(ConstrName,nv,1) repmat('  ',nv,1)...
            num2str((1:nv)','%d')  repmat('  ',nv,1) ...
            num2str([res.val1;res.val2;cstr.Min;cstr.Max;...
            res.penalty1;res.penalty1;]','%03.3e\t')]);
        ok=0;
    end

end
