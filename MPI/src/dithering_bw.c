#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "Util.h"

#define NB_BUFFERS 3

#define RIGHT 7 / 16
#define BOT_LEFT 3 / 16
#define BOT 5 / 16
#define BOT_RIGHT 1 / 16

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
    fprintf(output_file, "P5\n");
    fprintf(output_file, "%lu %lu\n", image->cols, image->rows);
    fprintf(output_file, "%lu\n", image->max_val);
    for (size_t i = 0; i < image->cols; i++) {
        for (size_t j = 0; j < image->rows; j++) {
            int16_t p = image->pixels[i * image->rows + j];
            if (p != 255 && p != 0)
                printf("ohoh p[i:%d, j:%d] = %d\n", j, i, p);
            // assert(p < 256 && p >= 0);
            // fprintf(output_file, "%d ", p);
            fprintf(output_file, "%c", p);
        }
        // fprintf(output_file, "\n");
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

void propagate_error_local(int16_t* data,
                           size_t rows,
                           size_t cols,
                           int16_t error,
                           size_t x,
                           size_t y) {
    if (x < cols - 1) {
        data[(y + 0) * cols + (x + 1)] =
            error * RIGHT + data[(y + 0) * cols + (x + 1)];
    }
    if (y < rows - 1) {
        if (x > 0) {
            data[(y + 1) * cols + (x - 1)] =
                error * BOT_LEFT + data[(y + 1) * cols + (x - 1)];
        }
        data[(y + 1) * cols + (x + 0)] =
            error * BOT + data[(y + 1) * cols + (x + 0)];
        if (x < cols - 1) {
            data[(y + 1) * cols + (x + 1)] =
                error * BOT_RIGHT + data[(y + 1) * cols + (x + 1)];
        }
    }
}

int16_t update_and_compute_error(int16_t* data, size_t i) {
    int16_t current_value = data[i];
    int16_t new_value = (current_value < 127) ? 0 : 255;
    // data[i] = new_value;
    *(data + i) = new_value;
    return current_value - new_value;
}

void sequential_floyd_steinberg(int16_t* data,
                                size_t rows,
                                size_t cols,
                                size_t offset_rows) {
    for (size_t y = 0; y < rows; y++) {
        for (size_t x = 0; x < cols; x++) {
            int16_t error = update_and_compute_error(data, y * cols + x);
            propagate_error_local(data, rows + offset_rows, cols, error, x, y);
        }
    }
}

void send_below(size_t block_index,
                size_t block_size,
                size_t line,
                size_t lines_per_process,
                int my_rank,
                int world_size,
                int16_t* error_to_bot) {
    if (block_index != 0 &&
        !(line == lines_per_process - 1 && my_rank == world_size - 1)) {
        MPI_Request req;
        // Error Buffer to send
        size_t index_block_to_send =
            ((block_index % NB_BUFFERS) + (NB_BUFFERS - 1)) % NB_BUFFERS;

        // printf("[%d] Sending to %d (line %d / %d) block #%d\n", my_rank,
        //        (my_rank + 1) % world_size, line, lines_per_process,
        //        block_index);
        size_t offset = index_block_to_send * block_size;
        MPI_Isend(error_to_bot + offset, block_size, MPI_INT16_T,
                  (my_rank + 1) % world_size, block_index - 1, MPI_COMM_WORLD,
                  &req);
        // TODO : ?
        // MPI_Request_free(&req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        // reinit values
        for (size_t i = 0; i < block_size; i++) {
            error_to_bot[i + offset] = 0;
        }
    }
}

void propagate_error_for_sending(int16_t* data,
                                 size_t block_index,
                                 size_t block_size,
                                 size_t blocks_per_line,
                                 size_t cols,
                                 size_t i,
                                 size_t offset,
                                 int16_t* error_to_bot,
                                 int16_t error) {
    // circular buffer index
    size_t base_index = (block_index % NB_BUFFERS) * block_size + i;
    size_t buf_size = NB_BUFFERS * block_size;

    if (block_index * block_size + i + 1 < cols) {
        data[offset + i + 1] = error * RIGHT + data[offset + i + 1];
        error_to_bot[(base_index + 1) % buf_size] += error * BOT_RIGHT;
    }
    error_to_bot[base_index % buf_size] += error * BOT;
    error_to_bot[(base_index - 1 + buf_size) % buf_size] += error * BOT_LEFT;
}

void fs_mpi_diagonal(int16_t* local_data,
                     size_t block_size,
                     size_t cols,
                     size_t lines_per_process,
                     size_t line_block_size) {
    int world_size = get_world_size();
    int my_rank = get_my_rank();

    /* ----- Local Dithering ----- */
    int16_t* error_from_top = calloc(block_size, sizeof(int16_t));
    size_t buf_size = NB_BUFFERS * block_size;
    int16_t* error_to_bot = calloc(buf_size, sizeof(int16_t));

    size_t y = 0;
    size_t blocks_per_line = cols / block_size;
    size_t iter = blocks_per_line + line_block_size - 1;
    assert(lines_per_process % line_block_size == 0);

    while (y < lines_per_process) {
        // printf("[%d] y = %lu\n", my_rank, y);
        for (size_t k = 0; k < iter; k++) {
            size_t i0 = (k >= blocks_per_line) ? blocks_per_line - 1 : k;
            size_t j0 = k - i0;
            // assert(i0 + j0 == k);
            /* ``i`` is in block_size */
            for (int32_t i = i0, j = j0; i >= 0 && j < line_block_size;
                 i--, j++) {
                size_t offset = (y + j) * cols + i * block_size;
                /* If I am not the top line of the image and if I am on the
                 * top line of the block of lines, I have to receive from
                 * the process above */
                if (j == 0 && !(my_rank == 0 && y == 0)) {
                    /* Receive the error from the process above */
                    MPI_Recv(error_from_top, block_size, MPI_INT16_T,
                             (my_rank + world_size - 1) % world_size, i,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    /* Update the pixels values with the errors received */
                    for (size_t b = 0; b < block_size; b++) {
                        local_data[offset + b] += error_from_top[b];
                    }
                }
                if (j == line_block_size - 1) {
                    /* We have to send to the process below */
                    for (size_t b = 0; b < block_size; b++) {
                        int16_t error =
                            update_and_compute_error(local_data, b + offset);

                        propagate_error_for_sending(
                            local_data, i, block_size, blocks_per_line, cols, b,
                            offset, error_to_bot, error);
                    }
                    send_below(i, block_size, y + j, lines_per_process, my_rank,
                               world_size, error_to_bot);
                } else {
                    for (size_t b = 0; b < block_size; b++) {
                        int16_t error =
                            update_and_compute_error(local_data, b + offset);

                        propagate_error_local(local_data, lines_per_process,
                                              cols, error, i * block_size + b,
                                              j + y);
                    }
                }
            }
        }
        y += line_block_size;
        send_below(blocks_per_line, block_size, y - 1, lines_per_process,
                   my_rank, world_size, error_to_bot);
    }
}

void fs_mpi(int16_t* local_data,
            size_t block_size,
            size_t cols,
            size_t lines_per_process,
            size_t line_block_size) {
    int world_size = get_world_size();
    int my_rank = get_my_rank();

    /* ----- Local Dithering ----- */
    int16_t* error_from_top = calloc(block_size, sizeof(int16_t));
    size_t buf_size = NB_BUFFERS * block_size;
    int16_t* error_to_bot = calloc(buf_size, sizeof(int16_t));
    size_t blocks_per_line = cols / block_size;
    // printf(
    //     "Process %d, lines_per_process: %d, line_block_size: %d, "
    //     "lines_per_process / line_block_size: %d, blocks_per_line %d \n",
    //     my_rank, lines_per_process, line_block_size,
    //     lines_per_process / line_block_size, blocks_per_line);

    // for (size_t block_of_lines = 0;
    //      block_of_lines < lines_per_process / line_block_size;
    //      block_of_lines++) {
    for (size_t line = 0; line < lines_per_process; line += line_block_size) {
        // printf("[%d] block_of_lines = %d\n", my_rank, block_of_lines);
        /* ----- First Line ----- */
        {
            for (size_t block_index = 0; block_index < blocks_per_line;
                 block_index++) {
                size_t block_offset = block_index * block_size;
                size_t offset = line * cols + block_offset;
                if (!(my_rank == 0 && line == 0)) {
                    /* ----- If this is not the top of the image, we receive the
                     * error from the above process */
                    // printf("[%d] Ready to recv from %d\n", my_rank,
                    //        (my_rank + world_size - 1) % world_size);
                    MPI_Recv(error_from_top, block_size, MPI_INT16_T,
                             (my_rank + world_size - 1) % world_size, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // printf("[%d] Done recving from %d\n", my_rank,
                    //        (my_rank + world_size - 1) % world_size);
                    /* ----- We add the error to the local data ----- */
                    for (size_t i = 0;
                         i < block_size && block_offset + i < cols; i++) {
                        local_data[offset + i] += error_from_top[i];
                    }
                }
                /* ----- Compute the local error ----- */
                for (size_t i = 0; i < block_size && block_offset + i < cols;
                     i++) {
                    int16_t error =
                        update_and_compute_error(local_data, i + offset);

                    if (line_block_size == 1) {
                        /* If we have to send right away */
                        propagate_error_for_sending(local_data, block_index,
                                                    block_size, blocks_per_line,
                                                    cols, i, offset,
                                                    error_to_bot, error);

                    } else {
                        propagate_error_local(
                            local_data, lines_per_process, cols, error,
                            i + block_size * block_index, line);
                    }
                    assert(local_data[i + offset] == 255 ||
                           local_data[i + offset] == 0);
                }
                /* If we have blocks of one line, we have to send the data to
                 * the process below */
                if (line_block_size == 1) {
                    send_below(block_index, block_size, line, lines_per_process,
                               my_rank, world_size, error_to_bot);
                }
            }
        }
        /* ----- Sequential Part ----- */
        if (line_block_size >= 3) {
            if (line_block_size > lines_per_process - line) {
                /* If the line_block_size is greater than the number of lines
                 * remaining, we update the value to make sure that we don't
                 * access junk data. */
                line_block_size = lines_per_process - line;
            }
            /* ``-2`` because we don't care about the first line of the block
             * (because we recv data) and the last line (because we send data).
             * So the only lines that we can do sequentially are between the
             * first and the last lines excluded */
            sequential_floyd_steinberg(local_data + (line + 1) * cols,
                                       line_block_size - 2, cols, 1);
            // printf("[%d] Done sequential part\n", my_rank);
        }
        // Last line
        if (line_block_size >= 2) {
            /* We just have to dither and send the data */
            size_t line_bot = line + line_block_size - 1;
            for (size_t block_index = 0; block_index < blocks_per_line;
                 block_index++) {
                size_t block_offset = block_index * block_size;
                size_t offset = line_bot * cols + block_offset;
                /* ----- Compute the local error ----- */
                for (size_t i = 0; i < block_size && block_offset + i < cols;
                     i++) {
                    int16_t error =
                        update_and_compute_error(local_data, i + offset);

                    propagate_error_for_sending(
                        local_data, block_index, block_size, blocks_per_line,
                        cols, i, offset, error_to_bot, error);

                    assert(local_data[i + offset] < 256 &&
                           local_data[i + offset] >= 0);
                }
                /* If we have blocks of one line, we have to send the data to
                 * the process below */
                send_below(block_index, block_size, line_bot, lines_per_process,
                           my_rank, world_size, error_to_bot);
            }
        }
        /* Sends the last block of the block of line */
        // printf("[%d] Ready to send the last block \n", my_rank);
        send_below(blocks_per_line, block_size, line + line_block_size - 1,
                   lines_per_process, my_rank, world_size, error_to_bot);
        // printf("[%d] Done sending the last block \n", my_rank);
    }
    // for (size_t i = 0; i < lines_per_process * cols; i++) {
    //     if (!(local_data[i] >= 0 && local_data[i] < 256)) {
    //         printf("[%d] local_data[%d] = %d ...\n", my_rank, i,
    //         local_data[i]);
    //     }
    // }
    free(error_from_top);
    free(error_to_bot);
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
            size_t block_offset = block_index * block_size;
            size_t offset = line_offset + block_offset;
            if (!(my_rank == 0 && line == 0)) {
                /* ----- If this is not the top of the image, we receive the
                 * error from the above process */
                MPI_Recv(error_from_top, block_size, MPI_INT16_T,
                         (my_rank + world_size - 1) % world_size, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            /* ----- We add the error to the local data ----- */
            for (size_t i = 0; i < block_size && block_offset + i < cols; i++) {
                local_data[offset + i] += error_from_top[i];
            }
            /* ----- Compute the local error ----- */
            for (size_t i = 0; i < block_size && block_offset + i < cols; i++) {
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
        for (size_t i = 0; i < block_size; i++) {
            // for (size_t i = 0; i < buf_size; i++) {
            (error_to_bot + ((blocks_per_line - 1) % 3) * block_size)[i] = 0;
            // error_to_bot[i] = 0;
        }
    }
    free(error_to_bot);
    free(error_from_top);
}

void floyd_steinberg(Image* image) {
    size_t rows = image->rows;
    size_t cols = image->cols;
    sequential_floyd_steinberg(image->pixels, rows, cols, 0);
}

size_t find_block_size(size_t w, size_t p) {
    size_t x = w / (2 * p);
    // printf("Block size: %ld\n", x);
    return x;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    char* filename = argv[1];
    char* out_filename = "out_mpi.pgm";
    uint32_t block_size = atoi(argv[2]);
    size_t line_block_size = atoi(argv[3]);

    Image* ppm_image;
    uint32_t h;
    uint32_t w;
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
        // block_size = (uint32_t)find_block_size(w, (size_t)world_size);
    }

    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    MPI_Bcast(&h, 1, MPI_UINT32_T, root, MPI_COMM_WORLD);
    MPI_Bcast(&w, 1, MPI_UINT32_T, root, MPI_COMM_WORLD);
    MPI_Bcast(&block_size, 1, MPI_UINT32_T, root, MPI_COMM_WORLD);

    /* ----- Computing the number of cells to send ----- */
    size_t lines_to_send_per_process = h / world_size;
    size_t cells_to_send_per_process = lines_to_send_per_process * w;
    int16_t* local_data = malloc(cells_to_send_per_process * sizeof(int16_t));

    MPI_Type_vector(lines_to_send_per_process / line_block_size,
                    w * line_block_size, world_size * w * line_block_size,
                    MPI_INT16_T, &PixelLine);
    MPI_Type_commit(&PixelLine);
    // printf("world_size %d, line_block_size %d\n", world_size,
    // line_block_size);
    /* ----- Sending the data ----- */
    if (my_rank == root) {
        MPI_Request req;
        for (size_t i = 0; i < world_size; i++) {
            MPI_Isend(pixels + i * w * line_block_size, 1, PixelLine, i, 0,
                      MPI_COMM_WORLD, &req);
        }
    }
    MPI_Recv(local_data, cells_to_send_per_process, MPI_INT16_T, root, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // floyd_steinberg_mpi(local_data, block_size, w,
    // lines_to_send_per_process);
    // fs_mpi(local_data, block_size, w, lines_to_send_per_process,
    //        line_block_size);
    fs_mpi_diagonal(local_data, block_size, w, lines_to_send_per_process,
                    line_block_size);
    // printf("[%d] Done dithering\n", my_rank);

    MPI_Request req;
    MPI_Isend(local_data, cells_to_send_per_process, MPI_INT16_T, root, 99,
              MPI_COMM_WORLD, &req);

    // printf("[%d] Done sending new pixels\n", my_rank);
    if (my_rank == root) {
        for (size_t i = 0; i < world_size; i++) {
            // printf("[%d] Ready to receive new pixels from %d\n", my_rank, i);
            MPI_Recv(pixels + i * w * line_block_size, 1, PixelLine, i, 99,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("[%d] Done receiving new pixels from %d\n", my_rank, i);
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

        write_image_to_file(ppm_image, "seq.pgm");

        free(ppm_image->pixels);
        free(ppm_image);
        double par_time = end_time - start_time;
        double seq_time = seq_end_time - seq_start_time;
        // H W P block_size r par seq
        printf("%d %d %d %d %lu %f %f\n", h, w, world_size, block_size,
               line_block_size, par_time, seq_time);
    } else {
        free(local_data);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
