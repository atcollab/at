classdef atoptions < handle
    %Definition of default parameters
    %
    % Singleton class for the storage of AT default values

    %% AT parameters
    properties
        XYStep = 3.0E-8             % Transverse step (for findorbitx, findmxx,...)
        DPStep = 3.0E-6             % Momentum step (for dispersion and chromaticity)
        OrbConvergence = 1.0E-12    % Convergence of findorbitx
        OrbMaxIter = 20             % Max. number of iterations for findorbitx
                                    % Number of OpenMP threads
        omp_num_threads = str2double(getenvopt('OMP_NUM_THREADS','0'))
    end
    
    %% AT dependent parameters
    properties (Dependent, SetAccess=private)
        openmp                      % True if compiled with OpenMP (not yet implemented)
    end
    
    methods
        function v=get.openmp(obj) %#ok<MANU>
            copts = coptions();
            v = logical(copts.openmp);
        end
    end
    
    %% options machinery
    methods (Static)
        function obj=instance(arg)
            % INSTANCE  Single instance of current AT options
            % INSTANCE()
            %   Return the single instance of current AT options
            %
            % INSTANCE('reset')
            %   Reset to the initial values and return the instance
            
            persistent vcurrent
            doreset=nargin > 0 && strcmp(arg,'reset');
            if doreset || isempty(vcurrent)
                vinitial=atoptions.initial(); %#ok<NASGU>
                vcurrent=atoptions();
            end
            obj=vcurrent;
        end
    end
    
    methods (Static, Access=private)
        function obj=initial()
            persistent vinitial
            if isempty(vinitial)
                vinitial=atoptions();
            end
            obj=vinitial;
        end
    end
    
    methods (Access=private)
        function obj=atoptions()
        end
    end
    
    methods
        function reset(obj, name)
            % RESET Revert the option to the initial value
            %   RESET(NAME) Reset the option NAME to the initial value
            obj.(name)=obj.initial().(name);
        end
    end
end
