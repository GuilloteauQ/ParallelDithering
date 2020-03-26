#include <omp.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "Util.h"

#define RIGHT 7 / 16
#define BOT_LEFT 3 / 16
#define BOT 5 / 16
#define BOT_RIGHT 1 / 16

#define CODE_RIGHT 8
#define CODE_BOT_LEFT 2
#define CODE_BOT 4
#define CODE_BOT_RIGHT 1

#define BLOCK_LINES 1

#define RANGE 256

#define BLOCK_SIZE 64

#define THRESHOLD 127

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
            fprintf(output_file, "%c", p);
        }
        // fprintf(output_file, "\n");
    }
    fclose(output_file);
}

Image* duplicate(Image* image, size_t iter) {
    size_t new_cols = image->cols * iter;
    size_t new_rows = image->rows * iter;

    int16_t* new_pixels = malloc(sizeof(int16_t) * new_cols * new_rows);

    for (size_t j = 0; j < new_rows; j++) {
        for (size_t i = 0; i < new_cols; i++) {
            new_pixels[i + new_cols * j] =
                image->pixels[(i % image->cols) +
                              (j % image->rows) * image->cols];
        }
    }

    Image* new_image = malloc(sizeof(Image));
    new_image->cols = new_cols;
    new_image->rows = new_rows;
    new_image->max_val = image->max_val;
    new_image->pixels = new_pixels;

    free(image->pixels);
    free(image);

    return new_image;
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
    // printf("data[%ld] = %d\n", i, *(data + i));
    int16_t current_value = *(data + i);
    int16_t new_value = (current_value < THRESHOLD) ? 0 : 255;
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

int16_t how_much_remaining(uint8_t remaining) {
    int16_t result = 0;
    for (uint8_t i = 1; i < 8; i += 2) {
        result += (remaining % 2 == 1) ? 8 * i : 0;
        remaining = remaining / 2;
    }
    return result;
}

int16_t random_value(uint8_t remaining) {
    int16_t result = 0;
    return 0;
    for (uint8_t i = 1; i < 8; i += 2) {
        result +=
            (remaining % 2 == 1) ? i * (rand() % RANGE - RANGE / 2) / 16 : 0;
        remaining = remaining / 2;
    }
    return result;
}

void init_tabs(uint8_t* remaining,
               pthread_mutex_t* locks,
               size_t rows,
               size_t cols) {
    for (size_t y = 0; y < rows; y++) {
        for (size_t x = 0; x < cols; x++) {
            size_t index = y * cols + x;
            pthread_mutex_init(&(locks[index]), NULL);
            if (x == 0)
                remaining[index] += CODE_RIGHT + CODE_BOT_RIGHT;
            if (x == cols - 1)
                remaining[index] += CODE_BOT_LEFT;
            if (y == 0)
                remaining[index] += CODE_BOT_LEFT + CODE_BOT + CODE_BOT_RIGHT;
        }
    }
}

void dither_block(int16_t* data,
                  size_t rows,
                  size_t cols,
                  uint8_t* remaining,
                  pthread_mutex_t* locks,
                  size_t x,
                  size_t y) {
    size_t X = x * BLOCK_SIZE;
    size_t reduced_cols = cols / BLOCK_SIZE;
    size_t index = y * cols + X;
    size_t max = ((x + 1) * BLOCK_SIZE < cols) ? BLOCK_SIZE
                                               : cols - (x + 0) * BLOCK_SIZE;
    pthread_mutex_t lock;
    for (size_t i = 0; i < max; i++) {
        int16_t error = update_and_compute_error(data, index + i);
        // Right: no need to lock because is in block
        if (X + i < cols - 1) {
            if (i == max - 1) {
                lock = locks[(y + 0) * reduced_cols + (x + 1)];
                // >>> LOCK
                pthread_mutex_lock(&lock);
                data[(y + 0) * cols + (X + i + 1)] += error * RIGHT;
                pthread_mutex_unlock(&lock);
                // <<< UNLOCK
            } else {
                data[(y + 0) * cols + (X + i + 1)] += error * RIGHT;
            }
        }
        if (y < rows - 1) {
            // Bot Left
            if (X + i > 0) {
                lock = (i == 0) ? locks[(y + 1) * reduced_cols + (x - 1)]
                                : locks[(y + 1) * reduced_cols + (x + 0)];
                // >>> LOCK
                pthread_mutex_lock(&lock);
                data[(y + 1) * cols + (X + i - 1)] += error * BOT_LEFT;
                pthread_mutex_unlock(&lock);
                // <<< UNLOCK
            }
            lock = locks[(y + 1) * reduced_cols + (x + 0)];
            // >>> LOCK
            pthread_mutex_lock(&lock);
            data[(y + 1) * cols + (X + i)] += error * BOT;
            pthread_mutex_unlock(&lock);
            // <<< UNLOCK
            if (X + i < cols - 1) {
                lock = (i == max - 1) ? locks[(y + 1) * reduced_cols + (x + 1)]
                                      : locks[(y + 1) * reduced_cols + (x + 0)];
                // >>> LOCK
                pthread_mutex_lock(&lock);
                data[(y + 1) * cols + (X + i + 1)] += error * BOT_RIGHT;
                pthread_mutex_unlock(&lock);
                // <<< UNLOCK
            }
        }
    }
    // Now updating the "remaining"
    if (x < reduced_cols - 1) {
        // We can update to our right
        index = (y + 0) * reduced_cols + (x + 1);
        lock = locks[index];
        // >>> LOCK
        pthread_mutex_lock(&lock);
        remaining[index] += CODE_RIGHT;
        if (remaining[index] >= 15) {
#pragma omp task
            dither_block(data, rows, cols, remaining, locks, x + 1, y + 0);
        }
        pthread_mutex_unlock(&lock);
        // <<< UNLOCK
    }
    if (y < rows - 1) {
        index = (y + 1) * reduced_cols + (x + 0);
        lock = locks[index];
        // >>> LOCK
        pthread_mutex_lock(&lock);
        remaining[index] += CODE_BOT;
        if (remaining[index] >= 15) {
#pragma omp task
            dither_block(data, rows, cols, remaining, locks, x + 0, y + 1);
        }
        pthread_mutex_unlock(&lock);
        // <<< UNLOCK
        if (x > 0) {
            index = (y + 1) * reduced_cols + (x - 1);
            lock = locks[index];
            // >>> LOCK
            pthread_mutex_lock(&lock);
            remaining[index] += CODE_BOT_LEFT;
            if (remaining[index] >= 15) {
#pragma omp task
                dither_block(data, rows, cols, remaining, locks, x - 1, y + 1);
            }
            pthread_mutex_unlock(&lock);
            // <<< UNLOCK
        }
        if (x < reduced_cols - 1) {
            index = (y + 1) * reduced_cols + (x + 1);
            lock = locks[index];
            // >>> LOCK
            pthread_mutex_lock(&lock);
            remaining[index] += CODE_BOT_RIGHT;
            if (remaining[index] >= 15) {
#pragma omp task
                dither_block(data, rows, cols, remaining, locks, x + 1, y + 1);
            }
            pthread_mutex_unlock(&lock);
            // <<< UNLOCK
        }
    }
#pragma omp taskwait
}

void parallel_floyd_steinberg_tasks_by_block(int16_t* data,
                                             size_t rows,
                                             size_t cols) {
    size_t n = rows * cols / BLOCK_SIZE;
    uint8_t* remaining = calloc(sizeof(uint8_t), n);
    pthread_mutex_t* locks = malloc(sizeof(pthread_mutex_t) * n);
    init_tabs(remaining, locks, rows, cols / BLOCK_SIZE);

    dither_block(data, rows, cols, remaining, locks, 0, 0);

    // for (size_t i = 0; i < n; i++) {
    //     pthread_mutex_destroy(&(locks[i]));
    // }
    // free(locks);
    // free(remaining);
}
void dither_pixel(int16_t* data,
                  size_t rows,
                  size_t cols,
                  uint8_t* remaining,
                  pthread_mutex_t* locks,
                  size_t x,
                  size_t y) {
    // printf("Dithering pixel (%ld, %ld) on thread %d\n", x, y,
    //        omp_get_thread_num());
    // I assume i have received every error
    size_t index = y * cols + x;
    int16_t error = update_and_compute_error(data, index);

    size_t index_error;
    size_t code_sum;
    if (y < rows - 1) {
        if (x > 0) {
            index_error = (y + 1) * cols + (x - 1);
            // >>> LOCK
            pthread_mutex_lock(&(locks[index_error]));
            data[index_error] += error * BOT_LEFT;
            remaining[index_error] += CODE_BOT_LEFT;
            code_sum = remaining[index_error];
            if (code_sum >= 15) {
#pragma omp task shared(data, remaining, locks, rows, cols, x, y)
                dither_pixel(data, rows, cols, remaining, locks, x - 1, y + 1);
            }
            pthread_mutex_unlock(&(locks[index_error]));
            // <<< UNLOCK
        }
        if (x < cols - 1) {
            index_error = (y + 1) * cols + (x + 1);
            // >>> LOCK
            pthread_mutex_lock(&(locks[index_error]));
            data[index_error] += error * BOT_RIGHT;
            remaining[index_error] += CODE_BOT_RIGHT;
            code_sum = remaining[index_error];
            if (code_sum >= 15) {
#pragma omp task shared(data, remaining, locks, rows, cols, x, y)
                dither_pixel(data, rows, cols, remaining, locks, x + 1, y + 1);
            }
            pthread_mutex_unlock(&(locks[index_error]));
            // <<< UNLOCK
        }
        index_error = (y + 1) * cols + (x + 0);
        // >>> LOCK
        pthread_mutex_lock(&(locks[index_error]));
        data[index_error] += error * BOT;
        remaining[index_error] += CODE_BOT;
        code_sum = remaining[index_error];
        if (code_sum >= 15) {
#pragma omp task shared(data, remaining, locks, rows, cols, x, y)
            dither_pixel(data, rows, cols, remaining, locks, x + 0, y + 1);
        }
        pthread_mutex_unlock(&(locks[index_error]));
        // <<< UNLOCK
    }
    if (x < cols - 1) {
        index_error = (y + 0) * cols + (x + 1);
        // >>> LOCK
        pthread_mutex_lock(&(locks[index_error]));
        data[index_error] += error * RIGHT;
        remaining[index_error] += CODE_RIGHT;
        code_sum = remaining[index_error];
        if (code_sum >= 15) {
#pragma omp task shared(data, remaining, locks, rows, cols, x, y)
            dither_pixel(data, rows, cols, remaining, locks, x + 1, y + 0);
        }
        pthread_mutex_unlock(&(locks[index_error]));
        // <<< UNLOCK
    }
}

void parallel_floyd_steinberg_tasks(int16_t* data, size_t rows, size_t cols) {
    size_t n = rows * cols;
    uint8_t* remaining = calloc(sizeof(uint8_t), n);
    pthread_mutex_t* locks = malloc(sizeof(pthread_mutex_t) * n);
    init_tabs(remaining, locks, rows, cols);

    dither_pixel(data, rows, cols, remaining, locks, 0, 0);

    // for (size_t i = 0; i < n; i++) {
    //     pthread_mutex_destroy(&(locks[i]));
    // }
    // free(locks);
    // free(remaining);
}

void parallel_floyd_steinberg(int16_t* data, size_t rows, size_t cols) {
    size_t n = rows * cols;
    uint8_t* remaining = calloc(sizeof(uint8_t), n);
    pthread_mutex_t* locks = malloc(sizeof(pthread_mutex_t) * n);
    init_tabs(remaining, locks, rows, cols);

    size_t y;
#pragma omp parallel for schedule(dynamic, BLOCK_LINES)
    for (y = 0; y < rows; y++) {
        for (size_t x = 0; x < cols; x++) {
            printf("y: %ld, x: %ld\n", y, x);
            size_t index = x + y * cols;
            pthread_mutex_t lock_current = locks[index];

            // >>> LOCK
            pthread_mutex_lock(&lock_current);
            int16_t current_value = data[index];
            int16_t error_max_remaining = how_much_remaining(remaining[index]);
            pthread_mutex_unlock(&lock_current);
            // <<< UNLOCK

            while (error_max_remaining < 15 &&
                   current_value + error_max_remaining >= THRESHOLD &&
                   current_value - error_max_remaining <= THRESHOLD) {
                printf("looping\n");
                // usleep(100);
                // >>> LOCK
                pthread_mutex_lock(&lock_current);
                current_value = data[index];
                error_max_remaining = how_much_remaining(remaining[index]);
                pthread_mutex_unlock(&lock_current);
                // <<< UNLOCK
            }
            // >>> LOCK
            pthread_mutex_lock(&lock_current);
            int16_t new_value = (current_value < THRESHOLD) ? 0 : 255;
            data[x + y * rows] = new_value;
            int16_t error =
                current_value + random_value(remaining[index]) - new_value;
            pthread_mutex_unlock(&lock_current);
            // <<< UNLOCK

            // printf("propagating error\n");
            size_t index_error;
            if (x < cols - 1) {
                index_error = (y + 0) * cols + (x + 1);
                // >>> LOCK
                pthread_mutex_lock(&(locks[index_error]));
                data[index_error] += error * RIGHT;
                remaining[index_error] += CODE_RIGHT;
                pthread_mutex_unlock(&(locks[index_error]));
                // <<< UNLOCK
            }
            if (y < rows - 1) {
                if (x > 0) {
                    index_error = (y + 1) * cols + (x - 1);
                    // >>> LOCK
                    pthread_mutex_lock(&(locks[index_error]));
                    data[index_error] += error * BOT_LEFT;
                    remaining[index_error] += CODE_BOT_LEFT;
                    pthread_mutex_unlock(&(locks[index_error]));
                    // <<< UNLOCK
                }
                index_error = (y + 1) * cols + (x + 0);
                // >>> LOCK
                pthread_mutex_lock(&(locks[index_error]));
                data[index_error] += error * BOT;
                remaining[index_error] += CODE_BOT;
                pthread_mutex_unlock(&(locks[index_error]));
                // <<< UNLOCK
                if (x < cols - 1) {
                    index_error = (y + 1) * cols + (x + 1);
                    // >>> LOCK
                    pthread_mutex_lock(&(locks[index_error]));
                    data[index_error] += error * BOT_RIGHT;
                    remaining[index_error] += CODE_BOT_RIGHT;
                    pthread_mutex_unlock(&(locks[index_error]));
                    // <<< UNLOCK
                }
            }
        }
    }
    free(remaining);
    for (size_t i = 0; i < n; i++) {
        pthread_mutex_destroy(&(locks[i]));
    }
    free(locks);
}

double floyd_steinberg_seq(Image* image) {
    size_t rows = image->rows;
    size_t cols = image->cols;
    double start, end;
#pragma omp parallel
    {
#pragma omp single
        start = omp_get_wtime();
        sequential_floyd_steinberg(image->pixels, rows, cols, 0);
        end = omp_get_wtime();
    }
    return end - start;
}

double floyd_steinberg_par(Image* image) {
    size_t rows = image->rows;
    size_t cols = image->cols;
    double start, end;
#pragma omp parallel
    {
#pragma omp single
        {
            start = omp_get_wtime();
            parallel_floyd_steinberg_tasks_by_block(image->pixels, rows, cols);
            end = omp_get_wtime();
        }
    }
    return end - start;
}
int main(int argc, char** argv) {
    srand(time(NULL));
    char* filename = argv[1];
    Image* ppm_image;
    ppm_image = read_image_from_file(filename);
    double par_time = floyd_steinberg_par(ppm_image);
    write_image_to_file(ppm_image, "par.pgm");

    ppm_image = read_image_from_file(filename);
    double seq_time = floyd_steinberg_seq(ppm_image);
    write_image_to_file(ppm_image, "seq.pgm");

    printf("Par time: %f, Seq time: %f, Speedup: %f\n", par_time, seq_time,
           seq_time / par_time);

    free(ppm_image->pixels);
    free(ppm_image);
    return 0;
}
