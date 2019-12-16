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

int main(int argc, char** argv) {
    char* filename = argv[1];
    char* out_filename = argv[2];
    size_t iter = atoi(argv[3]);
    Image* image = read_image_from_file(filename);
    image = duplicate(image, iter);
    write_image_to_file(image, out_filename);
    free(image->pixels);
    free(image);
    return 0;
}
