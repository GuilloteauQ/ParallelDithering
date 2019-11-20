#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "Util.h"

#define NB_BUFFERS 3
#define plop 0x4c000439

// taken from :
// https://stackoverflow.com/questions/40807833/sending-size-t-type-data-with-mpi
#if SIZE_MAX == UCHAR_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
#error "what is happening here?"
#endif

typedef struct image_t {
    size_t cols;
    size_t rows;
    size_t max_val;
    int16_t* pixels;
} Image;

Image* read_image_from_file(char* filename) {
    FILE* ppm_file = fopen(filename, "r");

    int ich1 = getc(ppm_file);
    if (ich1 == EOF)
        pm_erreur("EOF / read error / magic number");
    int ich2 = getc(ppm_file);
    int raw_byte = ich2 == '5';

    Image* image = malloc(sizeof(Image));

    image->cols = pm_getint(ppm_file);
    image->rows = pm_getint(ppm_file);
    image->max_val = pm_getint(ppm_file);

    image->pixels = malloc(sizeof(int16_t) * image->cols * image->rows);

    for (size_t i = 0; i < image->rows; i++) {
        for (size_t j = 0; j < image->cols; j++) {
            image->pixels[i * image->cols + j] = (int16_t)(
                (raw_byte) ? pm_getrawbyte(ppm_file) : pm_getint(ppm_file));
        }
    }

    fclose(ppm_file);
    return image;
}

void write_image_to_file(Image* image, char* filename) {
    FILE* output_file = fopen(filename, "w");
    fprintf(output_file, "P2\n");
    fprintf(output_file, "%lu %lu\n", image->cols, image->rows);
    fprintf(output_file, "%lu\n", image->max_val);
    for (size_t i = 0; i < image->cols; i++) {
        for (size_t j = 0; j < image->rows; j++) {
            int16_t p = image->pixels[i * image->rows + j];
            assert(p < 256 && p >= 0);
            fprintf(output_file, "%d ", p);
        }
        fprintf(output_file, "\n");
    }
    fclose(output_file);
}

int get_my_rank() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int get_world_size() {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

void floyd_steinberg_mpi(int16_t* local_data,
                         size_t block_size,
                         size_t cols,
                         size_t lines_per_process) {
    // assert(cols % block_size == 0);
    int world_size = get_world_size();
    int my_rank = get_my_rank();

    /* ----- Local Dithering ----- */
    int16_t* error_from_top = calloc(block_size, sizeof(int16_t));
    size_t buf_size = NB_BUFFERS * block_size;
    int16_t* error_to_bot = calloc(buf_size, sizeof(int16_t));
    size_t blocks_per_line = cols / block_size;

    for (size_t line = 0; line < lines_per_process; line++) {
        size_t line_offset = line * cols;
        for (size_t block_index = 0; block_index < blocks_per_line;
             block_index++) {
            size_t elem_offset = block_index * block_size;
            size_t offset = line_offset + elem_offset;
            if (!(my_rank == 0 && line == 0)) {
                /* ----- If this is not the top of the image, we receive the
                 * error from the above process */
                MPI_Recv(error_from_top, block_size, MPI_INT16_T,
                         (my_rank + world_size - 1) % world_size, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            /* ----- We add the error to the local data ----- */
            for (size_t i = 0; i < block_size && elem_offset + i < cols; i++) {
                local_data[offset + i] += error_from_top[i];
                // Should be useless ...
                // error_from_top[i] = 0;
            }
            /* ----- Compute the local error ----- */
            for (size_t i = 0; i < block_size && elem_offset + i < cols; i++) {
                int16_t current_value = local_data[offset + i];
                int16_t new_value = (current_value < 127) ? 0 : 255;
                local_data[offset + i] = new_value;
                int16_t error = current_value - new_value;
                // circular buffer index
                size_t base_index = (block_index % NB_BUFFERS) * block_size + i;

                if (!(block_index == blocks_per_line - 1 &&
                      i == block_size - 1)) {
                    local_data[offset + i + 1] =
                        error * 7 / 16 + local_data[offset + i + 1];

                    error_to_bot[(base_index + 1) % buf_size] += error * 1 / 16;
                }
                error_to_bot[base_index % buf_size] += error * 5 / 16;
                error_to_bot[(base_index - 1) % buf_size] += error * 3 / 16;
            }
            /* ----- We can send the previous errors ----- */
            if (block_index != 0 &&
                !(line == lines_per_process - 1 && my_rank == world_size - 1)) {
                MPI_Request req;
                size_t index_block_to_send =
                    ((block_index % NB_BUFFERS) + NB_BUFFERS - 1) % NB_BUFFERS;

                MPI_Isend(error_to_bot + index_block_to_send * block_size,
                          block_size, MPI_INT16_T, (my_rank + 1) % world_size,
                          0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                // reinit values
                for (size_t i = 0; i < block_size; i++) {
                    (error_to_bot + index_block_to_send * block_size)[i] = 0;
                }
            }
        }
        if (!(line == lines_per_process - 1 && my_rank == world_size - 1)) {
            MPI_Request req;
            MPI_Isend(error_to_bot +
                          ((blocks_per_line - 1) % NB_BUFFERS) * block_size,
                      block_size, MPI_INT16_T, (my_rank + 1) % world_size, 0,
                      MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
        }
        // reinit values
        // for (size_t i = 0; i < block_size; i++) {
        for (size_t i = 0; i < buf_size; i++) {
            // (error_to_bot + ((blocks_per_line - 1) % 3) * block_size)[i] = 0;
            error_to_bot[i] = 0;
        }
    }
    free(error_to_bot);
    free(error_from_top);
}

void floyd_steinberg(Image* image) {
    size_t rows = image->rows;
    size_t cols = image->cols;
    for (size_t y = 0; y < rows; y++) {
        for (size_t x = 0; x < cols; x++) {
            int16_t current_value = image->pixels[y * cols + x];
            int16_t new_value = (current_value < 127) ? 0 : 255;
            image->pixels[y * cols + x] = new_value;
            int16_t error = current_value - new_value;

            if (x < cols - 1) {
                image->pixels[(y + 0) * cols + (x + 1)] =
                    error * 7 / 16 + image->pixels[(y + 0) * cols + (x + 1)];
            }
            if (y < rows - 1) {
                if (x > 0) {
                    image->pixels[(y + 1) * cols + (x - 1)] =
                        error * 3 / 16 +
                        image->pixels[(y + 1) * cols + (x - 1)];
                }
                image->pixels[(y + 1) * cols + (x + 0)] =
                    error * 5 / 16 + image->pixels[(y + 1) * cols + (x + 0)];
                image->pixels[(y + 1) * cols + (x + 1)] =
                    error * 1 / 16 + image->pixels[(y + 1) * cols + (x + 1)];
            }
        }
    }
}

size_t find_block_size(size_t w, size_t p) {
    size_t x = w / (2 * p);
    printf("Block size: %ld\n", x);
    return x;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    char* filename = argv[1];
    char* out_filename = (argc == 3) ? argv[2] : "out_mpi.pgm";

    Image* ppm_image;
    uint32_t h;
    uint32_t w;
    uint32_t block_size;
    int world_size = get_world_size();
    int my_rank = get_my_rank();
    int root = world_size - 1;
    int16_t* pixels = NULL;
    MPI_Datatype PixelLine;

    if (my_rank == root) {
        ppm_image = read_image_from_file(filename);
        pixels = ppm_image->pixels;
        h = (uint32_t)ppm_image->rows;
        w = (uint32_t)ppm_image->cols;
        block_size = (uint32_t)find_block_size(w, (size_t)world_size);
    }
    // MPI_Bcast(&h, 1, MPI_UINT32_T, root, MPI_COMM_WORLD);
    // MPI_Bcast(&w, 1, MPI_UINT32_T, root, MPI_COMM_WORLD);
    // MPI_Bcast(&block_size, 1, MPI_UINT32_T, root, MPI_COMM_WORLD);

    MPI_Bcast(&h, 1, plop, root, MPI_COMM_WORLD);
    MPI_Bcast(&w, 1, plop, root, MPI_COMM_WORLD);
    MPI_Bcast(&block_size, 1, plop, root, MPI_COMM_WORLD);


    /* ----- Computing the number of cells to send ----- */
    size_t lines_to_send_per_process = h / world_size;
    size_t cells_to_send_per_process = lines_to_send_per_process * w;
    int16_t* local_data = malloc(cells_to_send_per_process * sizeof(int16_t));
    MPI_Type_vector(lines_to_send_per_process, w, world_size * w, MPI_INT16_T,
                    &PixelLine);
    MPI_Type_commit(&PixelLine);

    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    /* ----- Sending the data ----- */
    if (my_rank == root) {
        MPI_Request req;
        for (size_t i = 0; i < world_size; i++) {
            MPI_Isend(pixels + i * w, 1, PixelLine, i, 0, MPI_COMM_WORLD, &req);
        }
    }
    MPI_Recv(local_data, cells_to_send_per_process, MPI_INT16_T, root, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    floyd_steinberg_mpi(local_data, block_size, w, lines_to_send_per_process);

    MPI_Request req;
    MPI_Isend(local_data, cells_to_send_per_process, MPI_INT16_T, root, 0,
              MPI_COMM_WORLD, &req);

    if (my_rank == root) {
        for (size_t i = 0; i < world_size; i++) {
            MPI_Recv(pixels + i * w, 1, PixelLine, i, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }
    }
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    if (my_rank == root) {
        write_image_to_file(ppm_image, out_filename);

        double seq_start_time, seq_end_time;
        ppm_image = read_image_from_file(filename);
        seq_start_time = MPI_Wtime();
        floyd_steinberg(ppm_image);
        seq_end_time = MPI_Wtime();

        free(ppm_image->pixels);
        free(ppm_image);
        double par_time = end_time - start_time;
        double seq_time = seq_end_time - seq_start_time;
        double speedup = seq_time / par_time;
        double efficiency = speedup / (double)world_size;
        printf(
            "Execution Time:\n\tPar %f s\n\tSeq %f s\n\tSpeedup "
            "%f\n\tEfficiency %f\n",
            par_time, seq_time, speedup, efficiency);
    }
    free(local_data);
    MPI_Finalize();

    return 0;
}
