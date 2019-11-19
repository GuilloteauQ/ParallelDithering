#include <assert.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "Util.h"

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
    assert(cols % block_size == 0);
    int world_size = get_world_size();
    int my_rank = get_my_rank();

    /* ----- Local Dithering ----- */
    int16_t* error_from_top = calloc(block_size, sizeof(int16_t));
    int16_t* error_to_bot = calloc(block_size * 3, sizeof(int16_t));

    for (size_t line = 0; line < lines_per_process; line++) {
        for (size_t block_index = 0; block_index < cols / block_size;
             block_index++) {
            size_t offset = line * cols + block_index * block_size;
            if (!(my_rank == 0 && line == 0)) {
                /* ----- If this is not the top of the image, we receive the
                 * error from the above process */
                MPI_Recv(error_from_top, block_size, MPI_INT16_T,
                         (my_rank + world_size - 1) % world_size, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            /* ----- We add the error to the local data ----- */
            for (size_t i = 0; i < block_size; i++) {
                local_data[offset + i] += error_from_top[i];
                // Should be useless ...
                // error_from_top[i] = 0;
            }
            /* ----- Compute the local error ----- */
            for (size_t i = 0; i < block_size; i++) {
                int16_t current_value = local_data[offset + i];
                int16_t new_value = (current_value < 127) ? 0 : 255;
                local_data[offset + i] = new_value;
                int16_t error = current_value - new_value;
                // circular buffer index and size
                size_t base_index = (3 + (block_index % 3)) * block_size + i;
                size_t buf_size = 3 * block_size;

                if (!(block_index == cols / block_size &&
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
                MPI_Isend(error_to_bot +
                              (((block_index % 3) - 1 + 3) % 3) * block_size,
                          block_size, MPI_INT16_T, (my_rank + 1) % world_size,
                          0, MPI_COMM_WORLD, &req);
                MPI_Request_free(&req);
                // reinit values
                for (size_t i = 0; i < block_size; i++) {
                    (error_to_bot +
                     (((block_index % 3) - 1 + 3) % 3) * block_size)[i] = 0;
                }
            }
        }
        if (!(line == lines_per_process - 1 && my_rank == world_size - 1)) {
            MPI_Request req;
            MPI_Isend(error_to_bot + ((cols / block_size - 1) % 3) * block_size,
                      block_size, MPI_INT16_T, (my_rank + 1) % world_size, 0,
                      MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
        }
        // reinit values
        for (size_t i = 0; i < block_size; i++) {
            (error_to_bot + ((cols / block_size - 1) % 3) * block_size)[i] = 0;
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

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);


    char* filename = argv[1];
    size_t h = atoi(argv[2]);
    size_t w = atoi(argv[3]);
    size_t block_size = atoi(argv[4]);
    char* out_filename = (argc == 6) ? argv[5] : "out_mpi.pgm";

    Image* ppm_image;
    int world_size = get_world_size();
    int my_rank = get_my_rank();
    int16_t* pixels = NULL;
    MPI_Datatype PixelLine;

    /* ----- Computing the number of cells to send ----- */
    size_t lines_to_send_per_process = h / world_size;
    size_t cells_to_send_per_process = lines_to_send_per_process * w;
    int16_t* local_data = malloc(cells_to_send_per_process * sizeof(int16_t));
    MPI_Type_vector(lines_to_send_per_process, w, world_size * w, MPI_INT16_T,
                    &PixelLine);
    MPI_Type_commit(&PixelLine);

    if (my_rank == 0) {
        ppm_image = read_image_from_file(filename);
        pixels = ppm_image->pixels;
    }
    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    /* ----- Sending the data ----- */
    if (my_rank == 0) {
        MPI_Request req;
        for (size_t i = 0; i < world_size; i++) {
            MPI_Isend(pixels + i * w, 1, PixelLine, i, 0, MPI_COMM_WORLD, &req);
        }
    }
    MPI_Recv(local_data, cells_to_send_per_process, MPI_INT16_T, 0, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    floyd_steinberg_mpi(local_data, block_size, w, lines_to_send_per_process);

    MPI_Request req;
    MPI_Isend(local_data, cells_to_send_per_process, MPI_INT16_T, 0, 0,
              MPI_COMM_WORLD, &req);

    if (my_rank == 0) {
        for (size_t i = 0; i < world_size; i++) {
            MPI_Recv(pixels + i * w, 1, PixelLine, i, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }
    }
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    if (my_rank == 0) {
        write_image_to_file(ppm_image, out_filename);
        free(ppm_image->pixels);
        free(ppm_image);
	printf("Execution Time: %f s\n", end_time - start_time);
    }
    free(local_data);
    MPI_Finalize();

    return 0;
}
