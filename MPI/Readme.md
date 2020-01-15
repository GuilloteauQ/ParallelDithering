# Parallel Dithering: MPI Implementation

## How to build it ?

### Dependencies

You will need an implementation of MPI with ``mpicc`` and ``mpirun``.

### Build it

Just run: ``make``

## How to run it ?

Once the project has been built

```shell
mpirun -np 4 ./dithering_bw path/to/image.pgm output_filename
```

* 4: the number of processes

* ``path/to/image.pgm``: the path to the PGM image. It can be P2 or P5

* ``output_filename``: name of the output file (optional. default: ``out_mpi.pgm``)

## With Docker

### Download the Dockerfile

```bash
curl https://raw.githubusercontent.com/GuilloteauQ/ParallelDithering/master/MPI/Dockerfile -o Dockerfile
```

### Build the image

```bash
docker build -t mpi_dithering .
```

### Run

```bash
docker run mpi_dithering -p 16 -k 2 -r 64 -f 4
```

* ``-p``: number of processes (required)

* ``-k``: pixels per block (required)

* ``-r``: number of lines per block of lines (required)

* ``-f``: scaling factor of the image (integer)

* ``--image``: path to an image (default is mandrill)
