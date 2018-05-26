#include "ok_path.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int test_svg_parse() {
    ok_path_t *path1 = ok_path_alloc();
    ok_path_t *path2 = ok_path_alloc();

    // Test if two empty paths are equal
    if (!ok_path_equals(path1, path2)) {
        printf("Failure: %s: Empty paths not equal\n", __func__);
        ok_path_free(path1);
        ok_path_free(path2);
        return 1;
    }

    char *error;

    // Test if the SVG path is identical to the programatically constructed path.
    const char *svg_path = "M 100,100"
                           "L100,200 "
                           "h-100, " // Test superfluous comma
                           "v-100-100 " // Test no space or comma before negative number
                           "L0.25e-4,0.25E+2" // Test exponents
                           "M 0 0 a25,25 -30 0,1 50,-25 "
                           "M 200,300 Q400,50 600,300 T1000,300"
                           "M 100,200 C100,100 250,100 250,200 S400,300 400,200"
                           "Z";

    if (!ok_path_append_svg(path1, svg_path, &error)) {
        printf("Failure: %s: SVG parse error: %s\n", error, __func__);
        ok_path_free(path1);
        ok_path_free(path2);
        return 1;
    }

    ok_path_move_to(path2, 100, 100);
    ok_path_line_to(path2, 100, 200);
    ok_path_line_to(path2, 0, 200);
    ok_path_line_to(path2, 0, 100);
    ok_path_line_to(path2, 0, 0);
    ok_path_line_to(path2, 0.25e-4, 0.25e+2);
    ok_path_move_to(path2, 0, 0);
    ok_path_elliptical_arc_to(path2, 25, 25, -30 * M_PI / 180, false, true, 50, -25);
    ok_path_move_to(path2, 200, 300);
    ok_path_quad_curve_to(path2, 400, 50, 600, 300);
    ok_path_quad_curve_to(path2, 800, 550, 1000, 300);
    ok_path_move_to(path2, 100, 200);
    ok_path_curve_to(path2, 100, 100, 250, 100, 250, 200);
    ok_path_curve_to(path2, 250, 300, 400, 300, 400, 200);
    ok_path_close(path2);

    if (!ok_path_equals(path1, path2)) {
        printf("Failure: %s: Paths not equal\n", __func__);
        ok_path_free(path1);
        ok_path_free(path2);
        return 1;
    }

    printf("Success: %s\n", __func__);
    ok_path_free(path1);
    ok_path_free(path2);
    return 0;
}

static int test_iteration() {
    ok_path_t *path1 = ok_path_alloc();
    ok_path_t *path2 = ok_path_alloc();

    ok_path_move_to(path2, 100, 100);
    ok_path_line_to(path2, 100, 200);
    ok_path_line_to(path2, 0, 200);
    ok_path_line_to(path2, 0, 100);
    ok_path_line_to(path2, 0, 0);
    ok_path_line_to(path2, 0.25e-4, 0.25e+2);
    ok_path_move_to(path2, 0, 0);
    ok_path_elliptical_arc_to(path2, 25, 25, -30 * M_PI / 180, false, true, 50, -25);
    ok_path_move_to(path2, 200, 300);
    ok_path_quad_curve_to(path2, 400, 50, 600, 300);
    ok_path_quad_curve_to(path2, 800, 550, 1000, 300);
    ok_path_move_to(path2, 100, 200);
    ok_path_curve_to(path2, 100, 100, 250, 100, 250, 200);
    ok_path_curve_to(path2, 250, 300, 400, 300, 400, 200);
    ok_path_close(path2);

    ok_path_iterator_t iterator = 0;
    enum ok_path_segment_type type;
    double cx1, cy1;
    double cx2, cy2;
    double x, y;
    while (ok_path_segment_next(path2, &iterator, &type, &cx1, &cy1, &cx2, &cy2, &x, &y)) {
        switch (type) {
            case OK_PATH_MOVE_TO:
                ok_path_move_to(path1, x, y);
                break;
            case OK_PATH_LINE_TO:
                ok_path_line_to(path1, x, y);
                break;
            case OK_PATH_QUAD_CURVE_TO:
                ok_path_quad_curve_to(path1, cx1, cy1, x, y);
                break;
            case OK_PATH_CUBIC_CURVE_TO:
                ok_path_curve_to(path1, cx1, cy1, cx2, cy2, x, y);
                break;
            case OK_PATH_CLOSE:
                ok_path_close(path1);
                break;
        }
    }

    if (!ok_path_equals(path1, path2)) {
        printf("Failure: %s: Paths not equal\n", __func__);
        ok_path_free(path1);
        ok_path_free(path2);
        return 1;
    }

    printf("Success: %s\n", __func__);
    ok_path_free(path1);
    ok_path_free(path2);
    return 0;
}

static int test_append_lines() {
    ok_path_t *path1 = ok_path_alloc();
    ok_path_t *path2;

    char *error;
    const char *svg_path = "M 0,0 L 10,20, 40,20, 30,10, 50,0 Z";
    if (!ok_path_append_svg(path1, svg_path, &error)) {
        printf("Failure: %s: SVG parse error: %s\n", error, __func__);
        ok_path_free(path1);
        return 1;
    }

    // Using 2D array
    double points[][2] = {{0, 0}, {10, 20}, {40, 20}, {30, 10}, {50, 0}, {0, 0}};
    path2 = ok_path_alloc();
    ok_path_append_lines(path2, points, sizeof(points) / sizeof(*points));
    if (!ok_path_equals(path1, path2)) {
        printf("Failure: %s: Paths not equal\n", __func__);
        ok_path_free(path1);
        ok_path_free(path2);
        return 1;
    }

    // Using 1D array
    double points2[] = {0, 0, 10, 20, 40, 20, 30, 10, 50, 0, 0, 0};
    ok_path_reset(path2);
    ok_path_append_lines(path2, points, sizeof(points2) / (sizeof(*points2) * 2));
    if (!ok_path_equals(path1, path2)) {
        printf("Failure: %s: Paths not equal\n", __func__);
        ok_path_free(path1);
        ok_path_free(path2);
        return 1;
    }

    printf("Success: %s\n", __func__);
    ok_path_free(path1);
    ok_path_free(path2);
    return 0;
}

static int test_flatten() {
    ok_path_t *path = ok_path_alloc();

    char *error;

    const double start_x = 100;
    const double start_y = 100;
    const double end_x = 100;
    const double end_y = 200;
    const char *svg_path = "M 100,100 L100,200 h-100 v-100-100 L0.25e-4,0.25E+2"
        "M 0 0 a25,25 -30 0,1 50,-25 M 200,300 Q400,50 600,300 T1000,300 M 100,200 "
        "C100,100 250,100 250,200 S400,300 400,200 Z";

    if (!ok_path_append_svg(path, svg_path, &error)) {
        printf("Failure: %s: SVG parse error: %s\n", error, __func__);
        ok_path_free(path);
        return 1;
    }

    ok_flattened_path_t *flattened_path = ok_path_flatten(path);
    ok_path_free(path);

    double x, y;
    ok_flattened_path_location(flattened_path, 0.0, &x, &y, NULL);
    if (x != start_x || y != start_y) {
        printf("Failure: Flattened path error (start): %s\n", __func__);
        ok_flattened_path_free(flattened_path);
        return 1;
    }
    ok_flattened_path_location(flattened_path, 1.0, &x, &y, NULL);
    if (x != end_x || y != end_y) {
        printf("Failure: Flattened path error (end): %s\n", __func__);
        ok_flattened_path_free(flattened_path);
        return 1;
    }

    printf("Success: %s\n", __func__);
    ok_flattened_path_free(flattened_path);
    return 0;
}

int main() {
    return test_svg_parse() || test_iteration() || test_append_lines() || test_flatten();
}
