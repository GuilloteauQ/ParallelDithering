#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "Util.h"

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

int main(int argc, char** argv) {
    return 0;
}
