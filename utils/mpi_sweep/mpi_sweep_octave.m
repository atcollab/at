function mpi_sweep_octave(generator, worker, collector)

  pkg load mpi

  if not (MPI_Initialized ())
    MPI_Init ();
  endif

  comm = MPI_COMM_WORLD;
  info.rank = MPI_Comm_rank(comm);
  info.size = MPI_Comm_size(comm);
  info.name = MPI_Get_processor_name();

  N = 0;

  if info.rank == 0
    input = feval(generator);
    N = size(input)(2);
    array = split_array(input, info.size);
    for i=1:info.size
      MPI_Send(array{i}, i-1, 0, comm);
    endfor
  endif

  MPI_Barrier(comm);

  chunk = MPI_Recv(0, 0, comm);

  output = cell(1, size(chunk));
  for i=1:(size(chunk)(2))
    output{1, i} = {chunk{i}, worker(chunk{i}, info)};
  endfor

  MPI_Barrier(comm);

  MPI_Send(output, 0, 1, comm);

  MPI_Barrier(comm);

  if info.rank == 0
    array = cell(1, N);
    current = 1;
    for i=1:info.size
      chunk = MPI_Recv(i-1, 1, comm);
      for j=1:size(chunk)(2)
	array{current} = chunk{j};
	current = current + 1;
      endfor
    endfor
    collector(array)
  endif

  if not (MPI_Finalized)
    MPI_Finalize;
  endif

  function output=split_array(arr, num)
    indexes = zeros(num+1, 1);
    step = floor(size(arr)(2) / num);
    additional = size(arr)(2) - num * step;
    indexes(1) = 1;
    for i=2:num+1
      indexes(i) = indexes(i-1) + step;
      if additional > 0
	indexes(i) = indexes(i) + 1;
	additional = additional - 1;
      endif
    endfor
    output = cell(1, num);
    for i=1:num
      output{1,i} = arr(indexes(i):indexes(i+1)-1);
    endfor
  end
end
