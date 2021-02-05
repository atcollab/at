function varargout=test_openmp(threshold)
%TEST_OPENMP Test the perfoemance of OpenMP

ntot=50400;

% Compile without openMP
atmexall
[offa1, offa2] = run(ntot);
[offb1, offb2] = run(ntot);

% Compile with openMP, used foar any number of particles
atmexall -openmp -DOMP_PARTICLE_THRESHOLD=0
[ona1, ona2] = run(ntot);
[onb1, onb2] = run(ntot);
g1=(ona1+onb1)./(offa1+offb1);
g2=(ona2+onb2)./(offa2+offb2);

if nargin >= 1
    % Compile with openMP, used for more than 'threshold' particles
    arg=sprintf('-DOMP_PARTICLE_THRESHOLD=%i',threshold);
    atmexall('-openmp', arg);
    [opa1, opa2] = run(ntot);
    [opb1, opb2] = run(ntot);
    h1=(opa1+opb1)./(offa1+offb1);
    h2=(opa2+opb2)./(offa2+offb2);
    
    % Plot the result
    plot([g2 h2]);
    grid on
    title('Ratio with/without OpenMP');
    ylabel('T_{on}/T_{off}');
    xlabel('# particles');
    legend('{OMP\_PARTICLE\_THRESHOLD=0} ', sprintf('OMP\\_PARTICLE\\_THRESHOLD=%i',threshold));
    varargout={g1, g2, h1, h2};
else
    
    % Plot the result
    plot(g2);
    grid on
    title('Ratio with/without OpenMP');
    ylabel('T_{on}/T_{off}');
    xlabel('# particles');
    varargout={g1, g2};
end

    function [t1,t2]=run(nt)
        [t1,t2]=arrayfun(@(np) trackMultipleParticles('nparticles',np,'nturns',nt/np), (1:10)');
    end
end

