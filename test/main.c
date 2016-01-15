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
                           "h-100 "
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

    if (ok_path_get_length(path1) != ok_path_get_length(path2)) {
        printf("Failure: %s: Path lengths not equal\n", __func__);
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
    ok_path_free(path2);

    // Using 1D array
    double points2[] = {0, 0, 10, 20, 40, 20, 30, 10, 50, 0, 0, 0};
    path2 = ok_path_alloc();
    ok_path_append_lines(path2, points, sizeof(points2) / (sizeof(*points2) * 2));
    if (!ok_path_equals(path1, path2)) {
        printf("Failure: %s: Paths not equal\n", __func__);
        ok_path_free(path1);
        ok_path_free(path2);
        return 1;
    }
    ok_path_free(path2);

    printf("Success: %s\n", __func__);
    ok_path_free(path1);
    return 0;
}

int main() {
    return test_svg_parse() || test_append_lines();
}
