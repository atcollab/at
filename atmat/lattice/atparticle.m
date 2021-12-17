classdef atparticle
    %ATPARTICLE Particle definition for AT
    %   A particle is defined by its name, rest energy and charge

    properties(SetAccess=private)
        name
        rest_energy
        charge
    end

    properties(Constant, Hidden)
        e_mass = 1.0e6*PhysConstant.electron_mass_energy_equivalent_in_MeV.value;
        p_mass = 1.0e6*PhysConstant.proton_mass_energy_equivalent_in_MeV.value;
        known2 = struct(...
            'relativistic', {{0, -1}},...
            'electron', {{atparticle.e_mass, -1}},...
            'proton', {{atparticle.p_mass, 1}},...
            'positron', {{atparticle.e_mass, 1}}...
            );
    end

    methods(Static)
        function obj=loadobj(p)
            %LOADOBJ    Creates n ATPARTICLE from data in a .mat file
            obj = atparticle(p.name, p.rest_energy, p.charge);
        end
    end

    methods
        function obj = atparticle(name, varargin)
            %ATPARTICLE Construct an ATPARTICLE instance
            %   Detailed explanation goes here
            if ~strcmp(name,'electron')
                warning('AT:particle','AT tracking still assumes beta==1. Make sure your particle is ultra-relativistic');
            end
            try
                args = atparticle.known2.(name);
            catch
                args = varargin(1:2);
            end
            obj.name = name;
            [obj.rest_energy, obj.charge] = getargs(varargin, args{:});
        end

        function p = saveobj(obj)
            %SAVEOBJ    Generate the structure stored in a .mat file
            p = struct('name', obj.name, 'rest_energy', obj.rest_energy, 'charge', obj.charge);
        end
    end
end
