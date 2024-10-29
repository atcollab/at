classdef pytests < matlab.unittest.TestCase
    % ">> run(pytests)" to run all tests
    % ">> run(pytests, 'orbit4')" to run a the single orbit4 method

    properties(Constant)
        mlist=[...
            "machine_data/hmba",...
            "machine_data/dba",...
            "machine_data/spear3.m"];
    end
    
    properties(TestParameter)
        dp = {0., -0.01, 0.01};
        rad = struct("radoff","ring4","radon","ring6");
        lat = {"hmba","dba","spear3"};  % All lattices
        lat2 = {"hmba","spear3"};       % lattices with cavities
    end
    
    properties
        ring4
        ring6
    end

    methods(TestClassSetup)
        function load_lattice(testCase)
            % Shared setup for the entire test class
            t=warning('off','AT:atradon:NOCavity');
            setoption('WarningDp6D',false);
            for fpath=testCase.mlist
                [~,fname,~]=fileparts(fpath);
                [testCase.ring4.(fname),testCase.ring6.(fname)]=mload(fpath);
            end
            warning(t);

            function [ring4,ring6]=mload(fpath)
                mr=atradoff(atloadlattice(fullfile(atroot,'..',fpath)));
                pr=atwritepy(mr,'keep_all',true);
                ring4.m=mr;
                ring4.p=pr;
                ring6.m=atenable_6d(mr);
                ring6.p=pr.enable_6d(pyargs('copy',true));
            end
        end
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Static)
        function [dct,df]=dctdf(r4, dpp)
            [frf,l]=atGetRingProperties(r4,'rf_frequency','cell_length');
            [~,o0]=findorbit4(r4,dp=dpp,strict=-1);
            o1=ringpass(r4, o0);
            dct=o1(6);
            df=-frf*dct/(l+dct);
        end
    end

    methods(Test, TestTags="GitHub")
        % These tests may run in GitHub actions

        function lattice_pass(testCase,lat,rad)
            % Test of basic tracking
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
            [~,morbit4]=findorbit4(lattice.m,dp,'strict',-1);
            testCase.verifyEqual(morbit4,porbit4,AbsTol=1.E-15);
        end

        function syncorbit(testCase,lat,dp)
            lattice=testCase.ring4.(lat);
            [dct,~]=testCase.dctdf(lattice.m, dp);
            % python
            a=cell(lattice.p.find_sync_orbit(dct));
            [psyncorb,~]=deal(a{:});
            psyncorb=double(psyncorb)';
            % Matlab
            [~,msyncorb]=findsyncorbit(lattice.m,dct);
            testCase.verifyEqual(msyncorb,psyncorb,AbsTol=1.E-15);
        end

        function orbit6(testCase,lat2,dp)
            lattice=testCase.ring6.(lat2);
            % python
            a=cell(lattice.p.find_orbit6(pyargs(dp=dp)));
            [porbit6,~]=deal(a{:});
            porbit6=double(porbit6)';
            % Matlab
            [~,morbit6]=findorbit6(lattice.m,dp=dp);
            testCase.verifyEqual(morbit6,porbit6,AbsTol=2.E-12);
        end

        function m44(testCase,lat2,dp)
            lattice=testCase.ring4.(lat2);
            % Matlab
            mm44=findm44(lattice.m,dp);
            % python
            a2=cell(lattice.p.find_m44(dp));
            [pm44,~]=deal(a2{:});
            pm44=double(pm44);
            testCase.verifyEqual(mm44,pm44,AbsTol=3.e-9);
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

        function offmomdp(testCase, dp)
            % Checks that off-momentum is correctly taken into account
            % for 6D lattices
            [ring1, elem1]=atlinopt6(testCase.ring4.hmba.m,'get_chrom', dp=dp);
            [ring2, elem2]=atlinopt6(testCase.ring6.hmba.m,'get_chrom', dp=dp);
            testCase.verifyEqual(ring1.tune, ring2.tune(1:2), AbsTol=1e-6);
            testCase.verifyEqual(ring1.chromaticity, ring2.chromaticity(1:2), AbsTol=5.e-5);
            testCase.verifyEqual(elem1.beta, elem2.beta, AbsTol=1.e-4)
        end

        function offmomdct(testCase, dp)
            [dct,~]=testCase.dctdf(testCase.ring4.hmba.m, dp);
            [ring1, elem1]=atlinopt6(testCase.ring4.hmba.m,'get_chrom', dct=dct);
            [ring2, elem2]=atlinopt6(testCase.ring6.hmba.m,'get_chrom', dct=dct);
            testCase.verifyEqual(ring1.tune, ring2.tune(1:2), AbsTol=1e-6);
            testCase.verifyEqual(ring1.chromaticity, ring2.chromaticity(1:2), AbsTol=5.e-5);
            testCase.verifyEqual(elem1.beta, elem2.beta, AbsTol=1.e-4)
        end

        function offmomdf(testCase, dp)
            [~,df]=testCase.dctdf(testCase.ring4.hmba.m, dp);
            [ring1, elem1]=atlinopt6(testCase.ring4.hmba.m,'get_chrom', df=df);
            [ring2, elem2]=atlinopt6(testCase.ring6.hmba.m,'get_chrom', df=df);
            testCase.verifyEqual(ring1.tune, ring2.tune(1:2), AbsTol=1e-6);
            testCase.verifyEqual(ring1.chromaticity, ring2.chromaticity(1:2), AbsTol=5.e-5);
            testCase.verifyEqual(elem1.beta, elem2.beta, AbsTol=1.e-4)
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
            periodicity = atGetRingProperties(lattice.m,'Periodicity');
            [mtune,mchrom]=tunechrom(lattice.m,'get_chrom',dp=dp);
            ptune=double(lattice.p.get_tune(pyargs(dp=dp)));
            pchrom=double(lattice.p.get_chrom(pyargs(dp=dp)));
            testCase.verifyEqual(mod(mtune*periodicity,1),ptune,AbsTol=2.e-9);
            testCase.verifyEqual(mchrom*periodicity,pchrom,AbsTol=3.e-4);
        end

        function tunechrom6(testCase,lat2,dp)
            % test on and off-momentum tunes of 6D lattices
            lattice=testCase.ring6.(lat2);
            periodicity = atGetRingProperties(lattice.m,'Periodicity');
            mlat=atsetcavity(lattice.m,frequency='nominal',dp=dp);
            plat=lattice.p.set_rf_frequency(pyargs(dp=dp,copy=true));
            [mtune,mchrom]=tunechrom(mlat,'get_chrom');
            ptune=double(plat.get_tune());
            pchrom=double(plat.get_chrom());
            testCase.verifyEqual(mod(mtune*periodicity,1),ptune,AbsTol=1.e-9);
            testCase.verifyEqual(mchrom*periodicity,pchrom,RelTol=1.e-4,AbsTol=3.e-4);
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

        function avlin1(testCase, lat, dp)
            % Check average beta, dispersion, phase advance
            lattice=testCase.ring4.(lat);
            mrefs=true(1,length(lattice.m));
            prefs=py.numpy.array(mrefs);
            % python
            a=cell(lattice.p.avlinopt(dp, prefs));
            pbeta = double(a{2});
            pdisp = double(a{4});
            pmu = double(a{3});
            % Matlab
            [~,mbeta,mmu,mdisp,~,~]=atavedata(lattice.m,dp,mrefs);
            % check
            testCase.verifyEqual(mbeta,pbeta,AbsTol=1.E-8,RelTol=1.e-8);
            testCase.verifyEqual(mmu,pmu,AbsTol=1.E-8,RelTol=0);
            testCase.verifyEqual(mdisp,pdisp,AbsTol=1.E-8,RelTol=0);
        end

        function radiation_integrals(testCase,lat)
            lattice=testCase.ring4.(lat);
            mintegrals=atsummary(lattice.m,'NoDisplay').integrals(1:5);
            pintegrals=double(lattice.p.get_radiation_integrals());
            testCase.verifyEqual(mintegrals,pintegrals,RelTol=1.E-12);
        end

        function ringparameters(testCase,lat2)
            lattice=testCase.ring4.(lat2);
            [menergy,mharm,mperiodicity,mgamma,mmcf]=atGetRingProperties(lattice.m,...
                'Energy','HarmNumber','Periodicity','gamma','mcf');
            testCase.verifyEqual(menergy,lattice.p.energy);
            testCase.verifyEqual(mharm,double(lattice.p.harmonic_number));
            testCase.verifyEqual(mperiodicity,double(lattice.p.periodicity));
            testCase.verifyEqual(mgamma,lattice.p.gamma);
            testCase.verifyEqual(mmcf,lattice.p.mcf,RelTol=1.E-8);
        end

        function fastring(testCase,lat2)
            lattice=testCase.ring4.(lat2);
            [rm, rmrad]=atfastring(lattice.m);
            rpy=cell(py.at.fast_ring(lattice.p));
            rp=cell(rpy{1});
            rprad=cell(rpy{2});
            checkattr(rm, rp);
            checkattr(rmrad,rprad);

            function checkattr(rm, rp)
                testCase.verifyEqual(rm{2}.Frequency, rp{1}.Frequency, RelTol=1.0e-20);
                testCase.verifyEqual(rm{2}.Voltage, rp{1}.Voltage, RelTol=1.0e-20);
                testCase.verifyEqual(rm{3}.I2, rp{2}.I2, RelTol=1.0e-15);
                testCase.verifyEqual(rm{3}.Length, rp{2}.Length, RelTol=1.0e-20);
                testCase.verifyEqual(rm{3}.M66, double(rp{2}.M66), AbsTol=1.0e-7);
                testCase.verifyEqual(rm{end}.A1, rp{3}.A1, RelTol=0.01);
                testCase.verifyEqual(rm{end}.A2, rp{3}.A2, RelTol=0.02);
                testCase.verifyEqual(rm{end}.A3, rp{3}.A3, RelTol=0.01);
                testCase.verifyEqual(rm{end}.Alphax, rp{3}.Alphax, AbsTol=1.e-10);
                testCase.verifyEqual(rm{end}.Alphay, rp{3}.Alphay, AbsTol=1.e-10);
                testCase.verifyEqual(rm{end}.Betax, rp{3}.Betax, RelTol=1.e-10);
                testCase.verifyEqual(rm{end}.Betay, rp{3}.Betay, RelTol=1.e-10);
                testCase.verifyEqual(rm{end}.Qpx, rp{3}.Qpx, RelTol=1.e-8);
                testCase.verifyEqual(rm{end}.Qpy, rp{3}.Qpy, RelTol=1.e-8);
                if length(rm) >= 5
                    testCase.verifyEqual(rm{end-1}.Lmatp, double(rp{4}.Lmatp), AbsTol=2.e-7);
                end
            end
        end

        function emittances(testCase, lat2)
            % Check emittances, tunes and damping rates
            lattice=testCase.ring6.(lat2);
            % python
            pdata=cell(lattice.p.ohmi_envelope());
            pemit=double(py.getattr(pdata{2},'mode_emittances'));
            ptunes=double(py.getattr(pdata{2},'tunes'));
            pdamprate=double(py.getattr(pdata{2},'damping_rates'));
            %matlab
            [mdata,~,~,m]=ohmienvelope(lattice.m);
            jmt=jmat(3);
            aa=amat(m);
            nn=-aa'*jmt*mdata.R*jmt*aa;
            memit=0.5*[nn(1,1)+nn(2,2) nn(3,3)+nn(4,4) nn(5,5)+nn(6,6)];
            [mtunes,mdamprate]=atdampingrates(m);
            %check
            testCase.verifyEqual(memit,pemit,AbsTol=1.e-30,RelTol=1.e-6);
            testCase.verifyEqual(mtunes,ptunes,AbsTol=1.e-10);
            testCase.verifyEqual(mdamprate,pdamprate,AbsTol=1.e-10);
        end

    end
end
