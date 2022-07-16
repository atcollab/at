classdef pytests < matlab.unittest.TestCase

    properties(Constant)
        mlist=[...
            "pyat/machine_data/hmba",...
            "pyat/machine_data/dba",...
            "machine_data/spear3.m"];
    end
    
    properties
        ring4
        ring6
    end
    
    properties(TestParameter)
        dp = {0., -0.01, 0.01};
        dct = {0., -0.00005, 0.00005};
        lat = struct("hmba", "hmba","dba","dba","spear3","spear3");
        rad = struct("radoff","ring4","radon","ring6");
        lat2 = struct("hmba", "hmba","spear3","spear3");
    end

    methods(TestClassSetup)
        function load_lattice(testCase)
            % Shared setup for the entire test class
            t=warning('off','AT:atradon:NOCavity');
            for fpath=testCase.mlist
                [~,fname,~]=fileparts(fpath);
                [testCase.ring4.(fname),testCase.ring6.(fname)]=mload(fpath);
            end
            warning(t);

            function [ring4,ring6]=mload(fpath)
                mr=atradoff(atloadlattice(fullfile(atroot,'..',fpath)));
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

    methods(Test, TestTags="GitHub")
        % These tests may rub in GitHub actions

        function lattice_pass(testCase,lat,rad)
            % Test ob basic tracking
            lattice=testCase.(rad).(lat);
            rin=1.e-6*eye(6);
            pin=py.numpy.asfortranarray(rin);
            % python    Make a copy because the lattice_pass modifies its input
            pout=double(py.at.lattice_pass(lattice.p,pin.copy()));
            % Matlab
            mout=linepass(lattice.m,rin);
            testCase.verifyEqual(mout,pout,AbsTol=1.E-30);
        end

        function orbit4(testCase,lat,dp)
            lattice=testCase.ring4.(lat);
            % Python
            a=cell(lattice.p.find_orbit4(dp));
            [porbit4,~]=deal(a{:});
            porbit4=double(porbit4)';
            % Matlab
            [~,morbit4]=findorbit4(lattice.m,dp);
            testCase.verifyEqual(morbit4,porbit4,AbsTol=1.E-15);
        end

        function syncorbit(testCase,lat,dct)
            lattice=testCase.ring4.(lat);
            % python
            a=cell(lattice.p.find_sync_orbit(dct));
            [psyncorb,~]=deal(a{:});
            psyncorb=double(psyncorb)';
            % Matlab
            [~,msyncorb]=findsyncorbit(lattice.m,dct);
            testCase.verifyEqual(msyncorb,psyncorb,AbsTol=1.E-15);
        end

        function orbit6(testCase,lat2)
            lattice=testCase.ring6.(lat2);
            % python
            a=cell(lattice.p.find_orbit6());
            [porbit6,~]=deal(a{:});
            porbit6=double(porbit6)';
            % Matlab
            [~,morbit6]=findorbit6(lattice.m);
            testCase.verifyEqual(morbit6,porbit6,AbsTol=1.E-15 );
        end

        function m44(testCase,lat2,dp)
            lattice=testCase.ring4.(lat2);
            % Matlab
            mm44=findm44(lattice.m,dp);
            % python
            a2=cell(lattice.p.find_m44(dp));
            [pm44,~]=deal(a2{:});
            pm44=double(pm44);
            testCase.verifyEqual(mm44,pm44,AbsTol=2.e-9);
        end

        function m66(testCase,lat2)
            lattice=testCase.ring6.(lat2);
            % Matlab
            mm66=findm66(lattice.m);
            % python
            a2=cell(lattice.p.find_m66());
            [pm66,~]=deal(a2{:});
            pm66=double(pm66);
            testCase.verifyEqual(mm66,pm66,AbsTol=5.E-9);
        end
    end

    methods(Test)
        % These tests are disabled on GitHub because of a Lapack failure:
        % "Intel MKL ERROR: Parameter 11 was incorrect on entry to DGEEV."
        % preventing the computation if eigenvectors.
        % Hypothesis: library conflict when running python under Matlab

        function tunechrom4(testCase,lat,dp)
            % test on and off-momentum tunes of 4D lattices
            lattice=testCase.ring4.(lat);
            [mtune,mchrom]=tunechrom(lattice.m,'get_chrom',dp=dp);
            ptune=double(lattice.p.get_tune(pyargs(dp=dp)));
            pchrom=double(lattice.p.get_chrom(pyargs(dp=dp)));
            testCase.verifyEqual(mtune,ptune,AbsTol=2.e-9);
            testCase.verifyEqual(mchrom,pchrom,AbsTol=2.e-4);
        end

        function tunechrom6(testCase,lat2,dp)
            % test on and off-momentum tunes of 6D lattices
            lattice=testCase.ring6.(lat2);
            mlat=atsetcavity(lattice.m,frequency='nominal',dp=dp);
            plat=lattice.p.set_rf_frequency(pyargs(dp=dp,copy=true));
            [mtune,mchrom]=tunechrom(mlat,'get_chrom');
            ptune=double(plat.get_tune());
            pchrom=double(plat.get_chrom());
            testCase.verifyEqual(mtune,ptune,AbsTol=1.e-9);
            testCase.verifyEqual(mchrom,pchrom,AbsTol=1.e-4);
        end

        function linopt1(testCase,dp)
            %test linear optics of hmba 4D lattice
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
            testCase.verifyEqual(mtunes,ptunes,AbsTol=1.E-9);
            testCase.verifyEqual(mchrom,pchrom,AbsTol=2.E-5);
        end

        function radiation_integrals(testCase,lat)
            lattice=testCase.ring4.(lat);
            %mintegrals=atsummary(lattice.m,'NoDisplay').integrals(1:5);
            mintegrals=ringpara(lattice.m).integrals(1:5);
            pintegrals=double(lattice.p.get_radiation_integrals());
            testCase.verifyEqual(mintegrals,pintegrals,RelTol=1.E-12);
        end

        function ringparameters(testCase,lat2)
            lattice=testCase.ring4.(lat2);
            mprops=atGetRingProperties(lattice.m);
            mmcf=mcf(lattice.m,0.0);
            testCase.verifyEqual(mprops.Energy,lattice.p.energy);
            testCase.verifyEqual(mprops.HarmNumber,double(lattice.p.harmonic_number));
            testCase.verifyEqual(mprops.Periodicity,double(lattice.p.periodicity));
            testCase.verifyEqual(mmcf,lattice.p.mcf,RelTol=1.E-8);
        end
    end
end