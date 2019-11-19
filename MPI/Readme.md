# Parallel Dithering: MPI Implementation

## How to build it ?

### Dependencies

You will need an implementation of MPI with ``mpicc`` and ``mpirun``.

### Build it

Just run: ``make``

## How to run it ?

Once the project has been built

```shell
mpirun -np 4 ./dithering_bw path/to/image.pgm w h block_size output_filename
```

* 4: the number of processes

* ``path/to/image.pgm``: the path to the PGM image. It can be P2 or P5

* ``w``: width of the image

* ``h``: height of the image

* ``output_filename``: name of the output file (optional. default: ``out_mpi.pgm``)
