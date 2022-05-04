classdef pytests < matlab.unittest.TestCase

    properties
        ring4
        ring6
    end
    
    properties(TestParameter)
        dp = {0., -0.01, 0.01};
        dct = {0., 0.00005};
        lat = struct("hmba", "hmba","dba","dba");
        rad = struct("radoff","ring4","radon","ring6");
    end

    methods(TestClassSetup)
        function load_lattice(testCase)
            % Shared setup for the entire test class
            t=warning('off','AT:atradon:NOCavity');
            for fpath=["pyat/machine_data/hmba",...
                    "pyat/machine_data/dba"]
                [~,fname,~]=fileparts(fpath);
                [testCase.ring4.(fname),testCase.ring6.(fname)]=mload(fpath);
            end
            warning(t);

            function [ring4,ring6]=mload(fpath)
                mr=atloadlattice(fullfile(atroot,'..',fpath));
                pr=py.at.load_var(mr',pyargs('keep_all',true));
                ring4.m=mr;
                ring4.p=pr;
                ring6.m=atradon(mr);
                ring6.p=pr.radiation_on(pyargs('copy',true));
            end
        end
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        % Test methods

        function x_lattice_pass(testCase,lat,rad)
            lattice=testCase.(rad).(lat);
            rin=1.e-6*eye(6);
            pin=py.numpy.asfortranarray(rin);
            % python    Make a copy because the lattice_pass modifies its input
            pout=double(py.at.lattice_pass(lattice.p,pin.copy()));
            % Matlab
            mout=linepass(lattice.m,rin);
            testCase.verifyEqual(mout,pout,AbsTol=1.E-30);
        end

        function x_orbit4(testCase,lat,dp)
            lattice=testCase.ring4.(lat);
            % Python
            a=cell(lattice.p.find_orbit4(dp));
            [porbit4,~]=deal(a{:});
            porbit4=double(porbit4)';
            % Matlab
            [~,morbit4]=findorbit4(lattice.m,dp);
            testCase.verifyEqual(morbit4,porbit4,AbsTol=1.E-9);
        end

        function x_syncorbit(testCase,lat,dct)
            lattice=testCase.ring4.(lat);
            % python
            a=cell(lattice.p.find_sync_orbit(dct));
            [psyncorb,~]=deal(a{:});
            psyncorb=double(psyncorb)';
            % Matlab
            [~,msyncorb]=findsyncorbit(lattice.m,dct);
            testCase.verifyEqual(msyncorb,psyncorb,AbsTol=1.E-9);
        end

        function x_orbit6(testCase)
            lattice=testCase.ring6.hmba;
            % python
            a=cell(lattice.p.find_orbit6());
            [porbit6,~]=deal(a{:});
            porbit6=double(porbit6)';
            % Matlab
            [~,morbit6]=findorbit6(lattice.m);
            testCase.verifyEqual(morbit6,porbit6,AbsTol=1.E-9);
        end

        function linopt(testCase,dp)
            lattice=testCase.ring4.hmba;
            mrefs=true(1,length(lattice.m)+1);
            prefs=py.numpy.array(mrefs);
            % python
            a=cell(lattice.p.linopt6(prefs,pyargs(dp=dp,get_chrom=true)));
            [~,prdata,pedata]=deal(a{:});
            ptunes=double(py.getattr(prdata,'tune'));
            pchrom=double(py.getattr(prdata,'chromaticity'));
            pbeta=double(py.numpy.array(py.getattr(pedata,'beta')));
            % Matlab
            [mrdata,medata]=atlinopt6(lattice.m,mrefs,'get_chrom',dp=dp);
            mtunes=mrdata.tune;
            mchrom=mrdata.chromaticity;
            mbeta=cat(1,medata.beta);
            testCase.verifyEqual(mbeta,pbeta,AbsTol=1.E-9,RelTol=1.e-9);
            testCase.verifyEqual(mtunes,ptunes,AbsTol=1.E-9,RelTol=1.e-9);
            testCase.verifyEqual(mchrom,pchrom,AbsTol=2.E-5,RelTol=3.e-5);
        end
    end
end