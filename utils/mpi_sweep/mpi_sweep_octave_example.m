#!/usr/bin/env octave

function mpi_sweep_octave_example

  ## Make sure that AT source files are in path.
  ## For this go to atoctave folder and run
  ## > octave --eval 'bootstrap;savepath'

  ## this code will run on all MPI nodes
  D1.FamName = 'DR01';
  D1.Length  = 3;
  D1.PassMethod = 'DriftPass';

  QF.FamName = 'QF';
  QF.Length = 1;
  QF.K = 0.2;
  QF.PassMethod= 'QuadLinearPass';

  D2.FamName = 'DR02';
  D2.Length  = 3;
  D2.PassMethod = 'DriftPass';

  QD.FamName = 'QD';
  QD.Length = 1;
  QD.K = -0.2;
  QD.PassMethod= 'QuadLinearPass';

  FODOCELL = {D1 QF D2 QD};
  THERING = [FODOCELL];

  function output=generator()
    ## Generate input parameter list for worker function.
    ## It will only be executed on node 0 to save on allocations.
    output = {{0.1,-0.1}, {0.2,-0.2}, {0.3,-0.3}};
  end

  function output=worker(input, info)
    ## additional arguments may be used
    ## info.name;  # name of the processor
    ## info.rank;  # rank (id) of the node
    ## info.size;  # total number of nodes
    THERING{findcells(THERING,'FamName','QF')}.K = input{1,1};
    THERING{findcells(THERING,'FamName','QD')}.K = input{1,2};
    output = findm44(THERING,0);
  end

  function collector(input)
    ## Collect computed data.
    ## This function will only be executed on node 0.
    fid = fopen("output.csv", "w");
    fprintf(fid, "Kqf,Kqd");
    for i=1:4
      for j=1:4
        fprintf(fid, ",m%d%d", i, j);
      endfor
    endfor
    fprintf(fid, "\n");
    for ind=1:size(input)(2)
      params = input{1,ind}{1,1};
      output = input{1,ind}{1,2};
      fprintf(fid, "%d,%d", params{1,1}, params{1,2});
      for i=1:4
        for j=1:4
          fprintf(fid, ",%d", output(i,j));
        endfor
      endfor
      fprintf(fid, "\n");
    endfor
    fclose(fid);
  end

  mpi_sweep_octave(@generator, @worker, @collector);

end
