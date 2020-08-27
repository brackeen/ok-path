/*
 ok-path
 https://github.com/brackeen/ok-path
 Copyright (c) 2016-2020 David Brackeen

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "ok_path.h"
#include <math.h>
#include <memory.h>
#include <stdlib.h>

#ifndef NDEBUG
#include <stdio.h> // For snprintf
#endif

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

// MARK: Vector

// Examples:
// typedef struct vector_of(int) vector_int_t;
// struct vector_of_ints vector_of(int);
#define vector_of(type) \
    { type *values; size_t length; size_t capacity; }

#define vector_free(v) \
    free((v)->values)

#define vector_last(v) \
    ((v)->values + ((v)->length - 1))

#define vector_at(v, i) \
    ((v)->values + (i))

#define vector_push_new(v) \
    (vector_ensure_capacity(v, 1) ? ((v)->values + ((v)->length++)) : NULL)

#define vector_ensure_capacity(v, additional_count) \
    (((v)->length + (size_t)(additional_count) <= (v)->capacity) ? true : \
    vector_realloc((void **)&(v)->values, (v)->length + (size_t)(additional_count), \
                   sizeof(*(v)->values), &(v)->capacity))

static bool vector_realloc(void **values, size_t min_capacity, size_t element_size,
                           size_t *capacity) {
    size_t new_capacity = MAX(8, MAX(min_capacity, *capacity << 1));
    void *new_values = realloc(*values, element_size * new_capacity);
    if (new_values) {
        *values = new_values;
        *capacity = new_capacity;
        return true;
    } else {
        return false;
    }
}

// MARK: Path

static const double OK_PATH_DEFAULT_FLATNESS = 1.0;

#define ok_path_type_is_curve(type) ((type) == OK_PATH_QUAD_CURVE_TO || \
                                     (type) == OK_PATH_CUBIC_CURVE_TO)

struct ok_path_element {
    enum ok_path_element_type type;
    double x, y;

    // Control points. Valid if type is CURVE_TO, otherwise these are undefined.
    double cx1, cy1;
    double cx2, cy2;
};

struct ok_subpath {
    size_t first_index;
    size_t last_index;
    double origin_x;
    double origin_y;
    bool has_curves;
};

struct vector_of_path_elements vector_of(struct ok_path_element);

struct ok_path {
    struct vector_of_path_elements elements;
    struct vector_of(struct ok_subpath) subpaths;
    bool has_curves;
    double flatness;
};

static ok_path_t *ok_path_create_with_flatness(double flatness) {
    ok_path_t *path = calloc(1, sizeof(ok_path_t));
    path->flatness = flatness;
    return path;
}

ok_path_t *ok_path_create() {
    return ok_path_create_with_flatness(OK_PATH_DEFAULT_FLATNESS);
}

void ok_path_free(ok_path_t *path) {
    if (path) {
        vector_free(&path->elements);
        vector_free(&path->subpaths);
        free(path);
    }
}

double ok_path_get_flatness(const ok_path_t *path) {
    return path->flatness;
}

void ok_path_set_flatness(ok_path_t *path, double flatness) {
    path->flatness = MAX(0.01, flatness);
}

void ok_path_reset(ok_path_t *path) {
    path->elements.length = 0;
    path->subpaths.length = 0;
    path->has_curves = false;
}

static bool _ok_equals(double a, double b) {
    // Could be improved by using an epislon or something better.
    return (float)a == (float)b;
}

bool ok_path_equals(const ok_path_t *path1, const ok_path_t *path2) {

    if (path1->elements.length != path2->elements.length ||
        path1->subpaths.length != path2->subpaths.length) {
        return false;
    }
    struct ok_path_element *element1 = path1->elements.values;
    struct ok_path_element *element2 = path2->elements.values;
    for (size_t i = 0; i < path1->elements.length; i++) {
        enum ok_path_element_type type1 = element1->type;
        enum ok_path_element_type type2 = element2->type;
        if (type1 == OK_PATH_CLOSE) {
            type1 = OK_PATH_LINE_TO;
        }
        if (type2 == OK_PATH_CLOSE) {
            type2 = OK_PATH_LINE_TO;
        }

        if (type1 != type2) {
            return false;
        }

        switch (type1) {
            case OK_PATH_LINE_TO:
            case OK_PATH_MOVE_TO:
            case OK_PATH_CLOSE:
                if (!_ok_equals(element1->x, element2->x) ||
                    !_ok_equals(element1->y, element2->y)) {
                    return false;
                }
                break;
            case OK_PATH_QUAD_CURVE_TO:
                if (!_ok_equals(element1->cx1, element2->cx1) ||
                    !_ok_equals(element1->cy1, element2->cy1)) {
                    return false;
                }
                break;
            case OK_PATH_CUBIC_CURVE_TO:
                if (!_ok_equals(element1->cx1, element2->cx1) ||
                    !_ok_equals(element1->cy1, element2->cy1) ||
                    !_ok_equals(element1->cx2, element2->cx2) ||
                    !_ok_equals(element1->cy2, element2->cy2)) {
                    return false;
                }
                break;
        }

        element1++;
        element2++;
    }
    return true;
}

bool ok_path_is_flat(const ok_path_t *path) {
    return !path->has_curves;
}

size_t ok_path_element_count(const ok_path_t *path) {
    return path->elements.length;
}

enum ok_path_element_type ok_path_element_get(const ok_path_t *path, size_t index,
                                              double *out_cx1, double *out_cy1,
                                              double *out_cx2, double *out_cy2,
                                              double *out_x, double *out_y) {
    struct ok_path_element *element = path->elements.values + index;
    *out_x = element->x;
    *out_y = element->y;
    if (element->type == OK_PATH_QUAD_CURVE_TO) {
        *out_cx1 = element->cx1;
        *out_cy1 = element->cy1;
    } else if (element->type == OK_PATH_CUBIC_CURVE_TO) {
        *out_cx1 = element->cx1;
        *out_cy1 = element->cy1;
        *out_cx2 = element->cx2;
        *out_cy2 = element->cy2;
    }
    return element->type;
}

void ok_path_element_set(const ok_path_t *path, size_t index,
                         double cx1, double cy1, double cx2, double cy2, double x, double y) {
    struct ok_path_element *element = path->elements.values + index;
    element->cx1 = cx1;
    element->cy1 = cy1;
    element->cx2 = cx2;
    element->cy2 = cy2;
    element->x = x;
    element->y = y;
}

static double ok_path_last_x(const ok_path_t *path) {
    if (path->elements.length) {
        return vector_last(&path->elements)->x;
    } else {
        return 0;
    }
}

static double ok_path_last_y(const ok_path_t *path) {
    if (path->elements.length) {
        return vector_last(&path->elements)->y;
    } else {
        return 0;
    }
}

// MARK: Subpaths

size_t ok_subpath_count(const ok_path_t *path) {
    return path->subpaths.length;
}

ok_path_t *ok_subpath_create(const ok_path_t *path, size_t subpath_index) {
    ok_path_t *new_path = ok_path_create_with_flatness(path->flatness);
    const struct ok_subpath *subpath = vector_at(&path->subpaths, subpath_index);
    size_t count = subpath->last_index - subpath->first_index + 1;
    if (count > 0 && vector_ensure_capacity(&new_path->elements, count)) {
        struct ok_path_element *src_values = path->elements.values + subpath->first_index;
        struct ok_path_element *dst_values = new_path->elements.values;
        memcpy(dst_values, src_values, count * sizeof(struct ok_path_element));
        new_path->elements.length += count;
        new_path->has_curves = subpath->has_curves;

        struct ok_subpath *new_subpath = vector_push_new(&new_path->subpaths);
        new_subpath->first_index = 0;
        new_subpath->last_index = count - 1;
        new_subpath->origin_x = subpath->origin_x;
        new_subpath->origin_y = subpath->origin_y;
        new_subpath->has_curves = subpath->has_curves;
    } else {
        ok_path_free(new_path);
        new_path = NULL;
    }
    return new_path;
}

size_t ok_subpath_first_element_index(const ok_path_t *path, size_t subpath_index) {
    return vector_at(&path->subpaths, subpath_index)->first_index;
}

size_t ok_subpath_last_element_index(const ok_path_t *path, size_t subpath_index) {
    return vector_at(&path->subpaths, subpath_index)->last_index;
}

void ok_subpath_origin(const ok_path_t *path, size_t subpath_index, double *out_x, double *out_y) {
    *out_x = vector_at(&path->subpaths, subpath_index)->origin_x;
    *out_y = vector_at(&path->subpaths, subpath_index)->origin_y;
}

bool ok_subpath_is_flat(const ok_path_t *path, size_t subpath_index) {
    return !vector_at(&path->subpaths, subpath_index)->has_curves;
}

bool ok_subpath_is_closed(const ok_path_t *path, size_t subpath_index) {
    size_t element_index = ok_subpath_last_element_index(path, subpath_index);
    return vector_at(&path->elements, element_index)->type == OK_PATH_CLOSE;
}

// MARK: Modifying paths

static struct ok_path_element * _ok_path_new_element(ok_path_t *path,
                                                     enum ok_path_element_type type,
                                                     double x, double y) {
    // New subpath if current is MOVE, previous is CLOSE, or this is the first element
    if (type == OK_PATH_MOVE_TO) {
        struct ok_subpath *subpath = vector_push_new(&path->subpaths);
        subpath->first_index = path->elements.length;
        subpath->has_curves = false;
        subpath->origin_x = x;
        subpath->origin_y = y;
    } else if (path->elements.length == 0) {
        struct ok_subpath *subpath = vector_push_new(&path->subpaths);
        subpath->first_index = 0;
        subpath->has_curves = false;
        subpath->origin_x = 0.0;
        subpath->origin_y = 0.0;
    } else if (vector_last(&path->elements)->type == OK_PATH_CLOSE) {
        struct ok_subpath *subpath = vector_push_new(&path->subpaths);
        subpath->first_index = path->elements.length;
        subpath->has_curves = false;
        subpath->origin_x = vector_last(&path->elements)->x;
        subpath->origin_y = vector_last(&path->elements)->y;
    }
    struct ok_path_element *element = vector_push_new(&path->elements);
    element->type = type;
    element->x = x;
    element->y = y;
    vector_last(&path->subpaths)->last_index = path->elements.length - 1;
    if (ok_path_type_is_curve(type)) {
        path->has_curves = true;
        vector_last(&path->subpaths)->has_curves = true;
    }
    return element;
}

void ok_path_move_to(ok_path_t *path, double x, double y) {
    _ok_path_new_element(path, OK_PATH_MOVE_TO, x, y);
}

void ok_path_line_to(ok_path_t *path, double x, double y) {
    _ok_path_new_element(path, OK_PATH_LINE_TO, x, y);
}

void ok_path_quad_curve_to(ok_path_t *path, double cx, double cy, double x, double y) {
    struct ok_path_element *element = _ok_path_new_element(path, OK_PATH_QUAD_CURVE_TO, x, y);
    element->cx1 = cx;
    element->cy1 = cy;
}

void ok_path_curve_to(ok_path_t *path, double cx1, double cy1, double cx2, double cy2,
                      double x, double y) {
    struct ok_path_element *element = _ok_path_new_element(path, OK_PATH_CUBIC_CURVE_TO, x, y);
    element->cx1 = cx1;
    element->cy1 = cy1;
    element->cx2 = cx2;
    element->cy2 = cy2;
}

void ok_path_close(ok_path_t *path) {
    double x, y;
    if (path->subpaths.length == 0) {
        x = 0.0;
        y = 0.0;
    } else {
        x = vector_last(&path->subpaths)->origin_x;
        y = vector_last(&path->subpaths)->origin_y;
    }
    _ok_path_new_element(path, OK_PATH_CLOSE, x, y);
}

void ok_path_append(ok_path_t *path, const ok_path_t *path_to_append) {
    size_t count = path_to_append->elements.length;
    if (vector_ensure_capacity(&path->elements, count)) {
        // Append one element at a time, because elements are "commands" that may have different
        // meaning depending on the previous elements. For example, the first command of the
        // appended path could be line-to, which could (or could not) create a new subpath,
        // depending on context.
        for (size_t i = 0; i < count; i++) {
            struct ok_path_element *element = path_to_append->elements.values + i;
            switch (element->type) {
                case OK_PATH_MOVE_TO:
                    ok_path_move_to(path, element->x, element->y);
                    break;
                case OK_PATH_LINE_TO:
                    ok_path_line_to(path, element->x, element->y);
                    break;
                case OK_PATH_QUAD_CURVE_TO:
                    ok_path_quad_curve_to(path, element->cx1, element->cy1, element->x, element->y);
                    break;
                case OK_PATH_CUBIC_CURVE_TO:
                    ok_path_curve_to(path, element->cx1, element->cy1, element->cx2, element->cy2,
                                     element->x, element->y);
                    break;
                case OK_PATH_CLOSE:
                    ok_path_close(path);
                    break;
            }
        }
    }
}

void ok_path_append_lines(ok_path_t *path, enum ok_path_element_type first_type,
                          const double (*points)[2], size_t num_points) {
    struct vector_of_path_elements *path_elements = &path->elements;
    if (num_points > 0 && vector_ensure_capacity(path_elements, num_points)) {
        if (first_type == OK_PATH_MOVE_TO) {
            ok_path_move_to(path, (*points)[0], (*points)[1]);
        } else {
            ok_path_line_to(path, (*points)[0], (*points)[1]);
        }
        for (size_t i = 1; i < num_points; i++) {
            points++;
            ok_path_line_to(path, (*points)[0], (*points)[1]);
        }
    }
}

// MARK: Arc conversion

void ok_path_arc_to(ok_path_t *path, double radius, bool large_arc, bool sweep,
                    double x, double y) {
    ok_path_elliptical_arc_to(path, radius, radius, 0.0, large_arc, sweep, x, y);
}

void ok_path_elliptical_arc_to(ok_path_t *path, double radius_x, double radius_y,
                               double rotation_radians, bool large_arc,
                               bool sweep, double x, double y) {
    // Convert arc to a series of bezier curves.
    // Interface is the same as the SVG spec, except xAxisRotation is in radians, not degrees.
    // See http://www.w3.org/TR/SVG/paths.html#PathDataEllipticalArcCommands
    // See http://www.w3.org/TR/SVG/implnote.html#ArcImplementationNotes
    if (radius_x == 0.0 || radius_y == 0.0) {
        ok_path_line_to(path, x, y);
        return;
    }
    const double curr_x = ok_path_last_x(path);
    const double curr_y = ok_path_last_y(path);
    const double end_x = x;
    const double end_y = y;
    const double dx2 = (curr_x - end_x) / 2.0;
    const double dy2 = (curr_y - end_y) / 2.0;
    const double cos_a = cos(rotation_radians);
    const double sin_a = sin(rotation_radians);

    const double x1p = cos_a * dx2 + sin_a * dy2;
    const double y1p = -sin_a * dx2 + cos_a * dy2;
    const double x1p2 = x1p * x1p;
    const double y1p2 = y1p * y1p;

    // Correct out-of-range radii
    double rx = fabs(radius_x);
    double ry = fabs(radius_y);
    double rx2 = rx * rx;
    double ry2 = ry * ry;
    const double v = x1p2 / rx2 + y1p2 / ry2;
    if (v > 1.0) {
        const double v_sq = sqrt(v);
        rx = v_sq * rx;
        ry = v_sq * ry;
        rx2 = rx * rx;
        ry2 = ry * ry;
    }

    // Find center point (cx, cy)
    double S = sqrt(fmax(0, (rx2 * ry2 - rx2 * y1p2 - ry2 * x1p2) / (rx2 * y1p2 + ry2 * x1p2)));
    if (large_arc == sweep) {
        S = -S;
    }
    const double cxp = S * rx * y1p / ry;
    const double cyp = -S * ry * x1p / rx;
    const double cx = cos_a * cxp - sin_a * cyp + (curr_x + end_x) / 2.0;
    const double cy = sin_a * cxp + cos_a * cyp + (curr_y + end_y) / 2.0;

    // Find start_angle and end_angle
    const double ux = (x1p - cxp) / rx;
    const double uy = (y1p - cyp) / ry;
    const double vx = (-x1p - cxp) / rx;
    const double vy = (-y1p - cyp) / ry;
    double n = sqrt(ux * ux + uy * uy);
    double p = ux;
    double angle_start = acos(p / n);
    if (uy < 0.0) {
        angle_start = -angle_start;
    }
    n = sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy));
    p = ux * vx + uy * vy;
    double angle_extent = acos(p / n);
    if (ux * vy - uy * vx < 0.0) {
        angle_extent = -angle_extent;
    }
    if (!sweep && angle_extent > 0.0) {
        angle_extent -= 2.0 * M_PI;
    } else if (sweep && angle_extent < 0.0) {
        angle_extent += 2.0 * M_PI;
    }
    if (angle_extent == 0.0) {
        ok_path_line_to(path, x, y);
        return;
    }
    double angle_end = angle_start + angle_extent;

    // Create one bezier for each quadrant
    double cos_eta_b = cos(angle_start);
    double sin_eta_b = sin(angle_start);
    double a_cos_eta_b = rx * cos_eta_b;
    double b_sin_eta_b = ry * sin_eta_b;
    double a_sin_eta_b = rx * sin_eta_b;
    double b_cos_eta_b = ry * cos_eta_b;
    double x_b = cx + a_cos_eta_b * cos_a - b_sin_eta_b * sin_a;
    double y_b = cy + a_cos_eta_b * sin_a + b_sin_eta_b * cos_a;
    double x_b_dot = -a_sin_eta_b * cos_a - b_cos_eta_b * sin_a;
    double y_b_dot = -a_sin_eta_b * sin_a + b_cos_eta_b * cos_a;

    double prev_angle = angle_start;
    double d = (sweep ? M_PI_2 : -M_PI_2);
    bool done = false;
    int quad = 0;
    while (!done) {
        double angle = prev_angle + d;
        if ((++quad == 4) || (sweep && angle >= angle_end) || (!sweep && angle <= angle_end)) {
            angle = angle_end;
            done = true;
        }

        double da = angle - prev_angle;
        double alpha = 4.0 * tan(da / 4.0) / 3.0;
        prev_angle = angle;

        // Alternative alpha from: http://www.spaceroots.org/documents/ellipse/
        // XXX: test for very large arcs to see which alpha is better
        // double t = tan(da / 2.0);
        // double alpha = sin(da) * (fsqrt(4.0 + 3.0 * t * t) - 1.0) / 3.0;

        double x_a = x_b;
        double y_a = y_b;
        double x_a_dot = x_b_dot;
        double y_a_dot = y_b_dot;

        cos_eta_b = cos(angle);
        sin_eta_b = sin(angle);

        a_cos_eta_b = rx * cos_eta_b;
        b_sin_eta_b = ry * sin_eta_b;
        a_sin_eta_b = rx * sin_eta_b;
        b_cos_eta_b = ry * cos_eta_b;

        x_b_dot = -a_sin_eta_b * cos_a - b_cos_eta_b * sin_a;
        y_b_dot = -a_sin_eta_b * sin_a + b_cos_eta_b * cos_a;

        if (done) {
            x_b = end_x;
            y_b = end_y;
        } else {
            x_b = cx + a_cos_eta_b * cos_a - b_sin_eta_b * sin_a;
            y_b = cy + a_cos_eta_b * sin_a + b_sin_eta_b * cos_a;
        }

        ok_path_curve_to(path,
                         x_a + alpha * x_a_dot, y_a + alpha * y_a_dot,
                         x_b - alpha * x_b_dot, y_b - alpha * y_b_dot,
                         x_b, y_b);
    }
}

// MARK: SVG path parsing

static const char * const OK_PATH_SVG_COMMANDS = "MmLlHhVvAaQqTtCcSsZz";
static const unsigned int OK_PATH_SVG_COMMAND_VALUES[] = {
    2, 2, 2, 2, 1, 1, 1, 1, 7, 7,
    4, 4, 2, 2, 6, 6, 4, 4, 0, 0
};
static const char * const OK_PATH_SVG_WHITESPACE = " \t\r\n";
static const char * const OK_PATH_SVG_WHITESPACE_OR_COMMA = ", \t\r\n";

#ifdef NDEBUG
#define ok_path_error(out_error_message, format, ...) do { \
    if (out_error_message) *out_error_message = "ok_path_error"; \
} while (0)
#else
static char OK_PATH_SVG_ERROR_MESSAGE[80];
#define ok_path_error(out_error_message, format, ...) do { \
    if (out_error_message) { \
        snprintf(OK_PATH_SVG_ERROR_MESSAGE, sizeof(OK_PATH_SVG_ERROR_MESSAGE), \
                 format, __VA_ARGS__); \
        *out_error_message = OK_PATH_SVG_ERROR_MESSAGE; \
    } \
} while (0)
#endif

bool ok_path_append_svg(ok_path_t *path, const char *svg_path, char **out_error_message) {
    const char *str = svg_path;
    double values[7] = { 0 };
    char last_command = 0;
    unsigned int last_values_required = 0;
    bool last_control_set = false;
    double curr_x = 0.0, curr_y = 0.0;
    double control_x = 0.0, control_y = 0.0;

    while (str && *str) {
        bool command_required = last_values_required == 0;
        const char *skip_chars = (command_required ? OK_PATH_SVG_WHITESPACE :
                                  OK_PATH_SVG_WHITESPACE_OR_COMMA);
        if (strchr(skip_chars, *str)) {
            str++;
        } else {
            // Get command
            char command;
            unsigned int values_required;
            char *found = strchr(OK_PATH_SVG_COMMANDS, *str);
            if (found) {
                command = *str++;
                values_required = OK_PATH_SVG_COMMAND_VALUES[found - OK_PATH_SVG_COMMANDS];
            } else if (command_required) {
                command = *str++;
                values_required = 0;
            } else {
                command = last_command;
                values_required = last_values_required;
            }

            // Parse numbers
            if (values_required > 0) {
                unsigned int count = 0;
                while (*str && count < values_required) {
                    if (count > 0 && strchr(OK_PATH_SVG_WHITESPACE_OR_COMMA, *str)) {
                        str++;
                    } else {
                        char *endptr;
                        values[count++] = strtod(str, &endptr);
                        if (str == endptr) {
                            ok_path_error(out_error_message,
                                          "Could not parse number at: %li", (long)(str - svg_path));
                            return false;
                        }
                        str = endptr;
                    }
                }
                if (count < values_required) {
                    ok_path_error(out_error_message, "Unexpected EOF: Needed %i more numbers",
                                  (values_required - count));
                    return false;
                }
            }

            // Execute command
            bool control_set = false;
            switch (command) {
                case 'Z':
                case 'z': {
                    ok_path_close(path);
                    curr_x = vector_last(&path->subpaths)->origin_x;
                    curr_y = vector_last(&path->subpaths)->origin_y;
                    break;
                }
                case 'M': {
                    curr_x = values[0];
                    curr_y = values[1];
                    ok_path_move_to(path, curr_x, curr_y);
                    command = 'L';
                    break;
                }
                case 'm': {
                    curr_x += values[0];
                    curr_y += values[1];
                    ok_path_move_to(path, curr_x, curr_y);
                    command = 'l';
                    break;
                }
                case 'L': {
                    curr_x = values[0];
                    curr_y = values[1];
                    ok_path_line_to(path, curr_x, curr_y);
                    break;
                }
                case 'l': {
                    curr_x += values[0];
                    curr_y += values[1];
                    ok_path_line_to(path, curr_x, curr_y);
                    break;
                }
                case 'H': {
                    curr_x = values[0];
                    ok_path_line_to(path, curr_x, curr_y);
                    break;
                }
                case 'h': {
                    curr_x += values[0];
                    ok_path_line_to(path, curr_x, curr_y);
                    break;
                }
                case 'V': {
                    curr_y = values[0];
                    ok_path_line_to(path, curr_x, curr_y);
                    break;
                }
                case 'v': {
                    curr_y += values[0];
                    ok_path_line_to(path, curr_x, curr_y);
                    break;
                }
                case 'C': {
                    double c1x = values[0];
                    double c1y = values[1];
                    control_x = values[2];
                    control_y = values[3];
                    curr_x = values[4];
                    curr_y = values[5];
                    ok_path_curve_to(path, c1x, c1y, control_x, control_y, curr_x, curr_y);
                    control_set = true;
                    break;
                }
                case 'c': {
                    double c1x = curr_x + values[0];
                    double c1y = curr_y + values[1];
                    control_x = curr_x + values[2];
                    control_y = curr_y + values[3];
                    curr_x += values[4];
                    curr_y += values[5];
                    ok_path_curve_to(path, c1x, c1y, control_x, control_y, curr_x, curr_y);
                    control_set = true;
                    break;
                }
                case 'S':
                case 's': {
                    double c1x;
                    double c1y;
                    if (last_control_set) {
                        // "The reflection of the second control point on the previous command
                        // relative to the current point"
                        c1x = curr_x * 2.0 - control_x;
                        c1y = curr_y * 2.0 - control_y;
                    } else {
                        // "Coincident with the current point"
                        c1x = curr_x;
                        c1y = curr_y;
                    }
                    if (command == 'S') {
                        control_x = values[0];
                        control_y = values[1];
                        curr_x = values[2];
                        curr_y = values[3];
                    } else {
                        control_x = curr_x + values[0];
                        control_y = curr_y + values[1];
                        curr_x += values[2];
                        curr_y += values[3];
                    }
                    ok_path_curve_to(path, c1x, c1y, control_x, control_y, curr_x, curr_y);
                    control_set = true;
                    break;
                }
                case 'Q': {
                    control_x = values[0];
                    control_y = values[1];
                    curr_x = values[2];
                    curr_y = values[3];
                    ok_path_quad_curve_to(path, control_x, control_y, curr_x, curr_y);
                    control_set = true;
                    break;
                }
                case 'q': {
                    control_x = curr_x + values[0];
                    control_y = curr_y + values[1];
                    curr_x += values[2];
                    curr_y += values[3];
                    ok_path_quad_curve_to(path, control_x, control_y, curr_x, curr_y);
                    control_set = true;
                    break;
                }
                case 'T':
                case 't': {
                    if (last_control_set) {
                        // "The reflection of the second control point on the previous command
                        // relative to the current point"
                        control_x = curr_x * 2.0 - control_x;
                        control_y = curr_y * 2.0 - control_y;
                    } else {
                        // "Coincident with the current point"
                        control_x = curr_x;
                        control_y = curr_y;
                    }
                    if (command == 'T') {
                        curr_x = values[0];
                        curr_y = values[1];
                    } else {
                        curr_x += values[0];
                        curr_y += values[1];
                    }
                    ok_path_quad_curve_to(path, control_x, control_y, curr_x, curr_y);
                    control_set = true;
                    break;
                }
                case 'A':
                case 'a': {
                    double radius_x = values[0];
                    double radius_y = values[1];
                    double angle = values[2] * M_PI / 180.0;
                    bool large_arc = values[3] != 0.0;
                    bool sweep = values[4] != 0.0;
                    if (command == 'A') {
                        curr_x = values[5];
                        curr_y = values[6];
                    } else {
                        curr_x += values[5];
                        curr_y += values[6];
                    }
                    ok_path_elliptical_arc_to(path, radius_x, radius_y, angle, large_arc, sweep,
                                              curr_x, curr_y);
                    break;
                }
                default: {
                    ok_path_error(out_error_message, "Invalid SVG command %c at position %li",
                                  command, (long)(str - svg_path));
                    return false;
                }
            }
            last_command = command;
            last_values_required = values_required;
            last_control_set = control_set;
        }
    }
    return (str != NULL);
}

// MARK: Transforms

void ok_path_scale(ok_path_t *path, double scale_x, double scale_y) {
    for (size_t i = 0; i < ok_subpath_count(path); i++) {
        struct ok_subpath *subpath = path->subpaths.values + i;
        subpath->origin_x *= scale_x;
        subpath->origin_y *= scale_y;
    }
    for (size_t i = 0; i < ok_path_element_count(path); i++) {
        struct ok_path_element *element = path->elements.values + i;
        element->x *= scale_x;
        element->y *= scale_y;
        element->cx1 *= scale_x;
        element->cy1 *= scale_y;
        element->cx2 *= scale_x;
        element->cy2 *= scale_y;
    }
}

void ok_path_translate(ok_path_t *path, double translate_x, double translate_y) {
    for (size_t i = 0; i < ok_subpath_count(path); i++) {
        struct ok_subpath *subpath = path->subpaths.values + i;
        subpath->origin_x += translate_x;
        subpath->origin_y += translate_y;
    }
    for (size_t i = 0; i < ok_path_element_count(path); i++) {
        struct ok_path_element *element = path->elements.values + i;
        element->x += translate_x;
        element->y += translate_y;
        element->cx1 += translate_x;
        element->cy1 += translate_y;
        element->cx2 += translate_x;
        element->cy2 += translate_y;
    }
}

// MARK: Path flattening

typedef void (*ok_add_point_func)(enum ok_path_element_type type, double x, double y,
                                  void *userData);

static double _ok_approx_dist(double px, double py, double ax, double ay,
                              double bx, double by) {
    // Distance from a point to a line approximation. From Graphics Gems II, page 11.
    double dx = fabs(bx - ax);
    double dy = fabs(by - ay);

    double div = dx + dy - fmin(dx, dy) / 2.0;

    if (div == 0.0) {
        return 0.0;
    }

    double a2 = (py - ay) * (bx - ax) - (px - ax) * (by - ay);

    return fabs(a2) / div;
}

static size_t _ok_num_segments(double flatness,
                               double x0, double y0, double x1, double y1,
                               double x2, double y2, double x3, double y3) {
    double dist = fmax(_ok_approx_dist(x1, y1, x0, y0, x3, y3),
                       _ok_approx_dist(x2, y2, x0, y0, x3, y3));
    if (dist <= 0) {
        return 1;
    } else {
        // Found by trial and error when attempting to match Apple's Core Graphics.
        double num_segments = ceil(dist * (5 / (flatness + 1)));
        if (num_segments <= 1) {
            return 1;
        } else if (num_segments >= 256) {
            return 256; // XXX: Max of 256 segments? Why?
        } else {
            return (size_t)num_segments;
        }
    }
}

static void _ok_path_flatten_curve_division(ok_add_point_func add_point, void *userData,
                                            enum ok_path_element_type type,
                                            size_t num_segments,
                                            double x0, double y0, double x1, double y1,
                                            double x2, double y2, double x3, double y3) {
    const double t = 1.0 / num_segments;
    const double t2 = t * t;
    const double t3 = t2 * t;

    double xf = x0;
    double xfd = 3.0 * (x1 - x0) * t;
    double xfdd2 = 3.0 * (x0 - 2.0 * x1 + x2) * t2;
    double xfddd6 = (3.0 * (x1 - x2) + x3 - x0) * t3;
    double xfddd2 = 3.0 * xfddd6;
    double xfdd = 2.0 * xfdd2;
    double xfddd = 2.0 * xfddd2;

    double yf = y0;
    double yfd = 3.0 * (y1 - y0) * t;
    double yfdd2 = 3.0 * (y0 - 2.0 * y1 + y2) * t2;
    double yfddd6 = (3.0 * (y1 - y2) + y3 - y0) * t3;
    double yfddd2 = 3.0 * yfddd6;
    double yfdd = 2.0 * yfdd2;
    double yfddd = 2.0 * yfddd2;

    for (size_t i = 1; i < num_segments; i++) {
        xf += xfd + xfdd2 + xfddd6;
        xfd += xfdd + xfddd2;
        xfdd += xfddd;
        xfdd2 += xfddd2;

        yf += yfd + yfdd2 + yfddd6;
        yfd += yfdd + yfddd2;
        yfdd += yfddd;
        yfdd2 += yfddd2;

        add_point(type, xf, yf, userData);
    }

    add_point(type, x3, y3, userData);
}

static size_t _ok_path_flatten_curve(double flatness,
                                     ok_add_point_func add_point, void *userData,
                                     enum ok_path_element_type type,
                                     double x1, double y1, double x2, double y2,
                                     double x3, double y3, double x4, double y4) {
    // First division
    const double x12 = (x1 + x2) / 2;
    const double y12 = (y1 + y2) / 2;
    const double x23 = (x2 + x3) / 2;
    const double y23 = (y2 + y3) / 2;
    const double x34 = (x3 + x4) / 2;
    const double y34 = (y3 + y4) / 2;
    const double x123 = (x12 + x23) / 2;
    const double y123 = (y12 + y23) / 2;
    const double x234 = (x23 + x34) / 2;
    const double y234 = (y23 + y34) / 2;
    const double x1234 = (x123 + x234) / 2;
    const double y1234 = (y123 + y234) / 2;

    // Left division
    const double lx12 = (x1 + x12) / 2;
    const double ly12 = (y1 + y12) / 2;
    const double lx23 = (x12 + x123) / 2;
    const double ly23 = (y12 + y123) / 2;
    const double lx34 = (x123 + x1234) / 2;
    const double ly34 = (y123 + y1234) / 2;
    const double lx123 = (lx12 + lx23) / 2;
    const double ly123 = (ly12 + ly23) / 2;
    const double lx234 = (lx23 + lx34) / 2;
    const double ly234 = (ly23 + ly34) / 2;
    const double lx1234 = (lx123 + lx234) / 2;
    const double ly1234 = (ly123 + ly234) / 2;

    // Right division
    const double rx12 = (x1234 + x234) / 2;
    const double ry12 = (y1234 + y234) / 2;
    const double rx23 = (x234 + x34) / 2;
    const double ry23 = (y234 + y34) / 2;
    const double rx34 = (x34 + x4) / 2;
    const double ry34 = (y34 + y4) / 2;
    const double rx123 = (rx12 + rx23) / 2;
    const double ry123 = (ry12 + ry23) / 2;
    const double rx234 = (rx23 + rx34) / 2;
    const double ry234 = (ry23 + ry34) / 2;
    const double rx1234 = (rx123 + rx234) / 2;
    const double ry1234 = (ry123 + ry234) / 2;

    // Determine the number of segments for each division
    size_t num_segments1 = _ok_num_segments(flatness, x1, y1, lx12, ly12, lx123, ly123, lx1234, ly1234);
    size_t num_segments2 = _ok_num_segments(flatness, lx1234, ly1234, lx234, ly234, lx34, ly34, x1234, y1234);
    size_t num_segments3 = _ok_num_segments(flatness, x1234, y1234, rx12, ry12, rx123, ry123, rx1234, ry1234);
    size_t num_segments4 = _ok_num_segments(flatness, rx1234, ry1234, rx234, ry234, rx34, ry34, x4, y4);

    // Convert to lines
    if (add_point) {
        _ok_path_flatten_curve_division(add_point, userData, type, num_segments1,
                                        x1, y1, lx12, ly12, lx123, ly123, lx1234, ly1234);
        _ok_path_flatten_curve_division(add_point, userData, type, num_segments2,
                                        lx1234, ly1234, lx234, ly234, lx34, ly34, x1234, y1234);
        _ok_path_flatten_curve_division(add_point, userData, type, num_segments3,
                                        x1234, y1234, rx12, ry12, rx123, ry123, rx1234, ry1234);
        _ok_path_flatten_curve_division(add_point, userData, type, num_segments4,
                                        rx1234, ry1234, rx234, ry234, rx34, ry34, x4, y4);
    }
    return num_segments1 + num_segments2 + num_segments3 + num_segments4;
}

static size_t _ok_path_flatten_generic(const ok_path_t *path,
                                       size_t first_subpath, size_t last_subpath,
                                       bool normalize_subpaths, bool close_subpaths,
                                       ok_add_point_func add_point, void *userData) {
    size_t count = 0;
    for (size_t subpath_index = first_subpath; subpath_index <= last_subpath; subpath_index++) {
        size_t first_index = ok_subpath_first_element_index(path, subpath_index);
        size_t last_index = ok_subpath_last_element_index(path, subpath_index);
        
        double x, y;
        ok_subpath_origin(path, subpath_index, &x, &y);

        if (normalize_subpaths && (path->elements.values + first_index)->type != OK_PATH_MOVE_TO) {
            if (add_point) {
                add_point(OK_PATH_MOVE_TO, x, y, userData);
            }
            count++;
        }

        if (ok_subpath_is_flat(path, subpath_index)) {
            if (add_point) {
                for (size_t i = first_index; i <= last_index; i++) {
                    struct ok_path_element *element = &path->elements.values[i];
                    add_point(element->type, element->x, element->y, userData);
                }
            }
            count += last_index - first_index + 1;
        } else {
            for (size_t i = first_index; i <= last_index; i++) {
                struct ok_path_element *element = &path->elements.values[i];
                if (element->type == OK_PATH_QUAD_CURVE_TO) {
                    double cx1 = (x + element->cx1 * 2.0) / 3.0;
                    double cy1 = (y + element->cy1 * 2.0) / 3.0;
                    double cx2 = (element->x + element->cx1 * 2.0) / 3.0;
                    double cy2 = (element->y + element->cy1 * 2.0) / 3.0;
                    count += _ok_path_flatten_curve(path->flatness, add_point, userData,
                                                    OK_PATH_QUAD_CURVE_TO, x, y,
                                                    cx1, cy1, cx2, cy2, element->x, element->y);
                } else if (element->type == OK_PATH_CUBIC_CURVE_TO) {
                    count += _ok_path_flatten_curve(path->flatness, add_point, userData,
                                                    OK_PATH_CUBIC_CURVE_TO, x, y,
                                                    element->cx1, element->cy1,
                                                    element->cx2, element->cy2,
                                                    element->x, element->y);
                } else {
                    if (add_point) {
                        add_point(element->type, element->x, element->y, userData);
                    }
                    count++;
                }
                x = element->x;
                y = element->y;
            }
        }
        if (close_subpaths && !ok_subpath_is_closed(path, subpath_index)) {
            if (add_point) {
                struct ok_path_element *element = &path->elements.values[first_index];
                add_point(OK_PATH_CLOSE, element->x, element->y, userData);
            }
            count++;
        }
    }
    return count;
}

static void _ok_path_add_point(enum ok_path_element_type type, double x, double y, void *userData) {
    ok_path_t *flattened_path = userData;
    if (type == OK_PATH_MOVE_TO) {
        ok_path_move_to(flattened_path, x, y);
    } else if (type == OK_PATH_CLOSE) {
        ok_path_close(flattened_path);
    } else {
        ok_path_line_to(flattened_path, x, y);
    }
}

static ok_path_t *_ok_path_flatten(const ok_path_t *path, size_t first_subpath,
                                   size_t last_subpath, bool close_subpaths) {
    size_t count = _ok_path_flatten_generic(path, first_subpath, last_subpath, false,
                                            close_subpaths, NULL, NULL);
    ok_path_t *flattened_path = ok_path_create_with_flatness(path->flatness);
    if (!vector_ensure_capacity(&flattened_path->elements, count)) {
        ok_path_free(flattened_path);
        return NULL;
    } else {
        _ok_path_flatten_generic(path, first_subpath, last_subpath, false, close_subpaths,
                                 _ok_path_add_point, flattened_path);
        return flattened_path;
    }
}

ok_path_t *ok_path_flatten(const ok_path_t *path) {
    if (path->has_curves) {
        return _ok_path_flatten(path, 0, path->subpaths.length - 1, false);
    } else {
        ok_path_t *new_path = ok_path_create_with_flatness(path->flatness);
        ok_path_append(new_path, path);
        return new_path;
    }
}

ok_path_t *ok_subpath_flatten(const ok_path_t *path, size_t subpath_index) {
    const struct ok_subpath *subpath = vector_at(&path->subpaths, subpath_index);
    if (subpath->has_curves) {
        return _ok_path_flatten(path, subpath_index, subpath_index, false);
    } else {
        return ok_subpath_create(path, subpath_index);
    }
}

struct ok_point_list_context {
    void *buffer;
    size_t offset;
    size_t stride;
    size_t count;
};

static void _ok_subpath_add_point_to_list(enum ok_path_element_type type, double x, double y,
                                          void *userData) {
    (void)type;
    struct ok_point_list_context *context = userData;
    double point[2] = { x, y };
    size_t offset = context->offset + context->stride * context->count;
    if (context->count > 0) {
        // Don't add duplicate points
        void *prev_point = (uint8_t *)context->buffer + (offset - context->stride);
        if (memcmp(point, prev_point, sizeof(point)) == 0) {
            return;
        }
    }
    void *dst = (uint8_t *)context->buffer + offset;
    memcpy(dst, point, sizeof(point));
    context->count++;
}

void ok_subpath_create_point_list_generic(const ok_path_t *path, size_t subpath_index,
                                          size_t offset, size_t stride,
                                          void **out_points, size_t *out_num_points) {
    if (stride < offset + sizeof(double[2])) {
        *out_points = NULL;
        *out_num_points = 0;
        return;
    }
    size_t count = _ok_path_flatten_generic(path, subpath_index, subpath_index, true, false,
                                            NULL, NULL);
    void *buffer = malloc(count * stride);
    if (!buffer) {
        *out_points = NULL;
        *out_num_points = count;
        return;
    }
    struct ok_point_list_context context = {
        .buffer = buffer, .offset = offset, .stride = stride, .count = 0 };
    _ok_path_flatten_generic(path, subpath_index, subpath_index, true, false,
                             _ok_subpath_add_point_to_list, &context);
    
    // NOTE: This would be a waste of memory if the path has many consecutive duplicate points.
    // However, this is an unlikely scenario.
    *out_points = buffer;
    *out_num_points = context.count;
}

// MARK: Motion Paths

struct ok_motion_path_segment {
    // The type of command this line segment originated from.
    enum ok_path_element_type type;

    // Point location.
    double x, y;

    // Entire length of the path up to this point.
    double length_to;

    // Angle from previous point to this point.
    double angle_to;
};

struct ok_motion_path {
    struct vector_of(struct ok_motion_path_segment) segments;
};

void ok_motion_path_free(ok_motion_path_t *path) {
    if (path) {
        vector_free(&path->segments);
        free(path);
    }
}

static void _ok_motion_path_add_point(enum ok_path_element_type type, double x, double y,
                                      void *userData) {
    ok_motion_path_t *path = userData;
    double dx;
    double dy;
    double prev_length;
    if (path->segments.length == 0) {
        dx = x;
        dy = y;
        prev_length = 0.0;
    } else {
        struct ok_motion_path_segment *segment = vector_last(&path->segments);
        dx = x - segment->x;
        dy = y - segment->y;
        prev_length = segment->length_to;
    }

    struct ok_motion_path_segment *segment = vector_push_new(&path->segments);
    if (segment) {
        segment->type = type;
        segment->x = x;
        segment->y = y;
        segment->angle_to = atan2(dy, dx);
        if (type == OK_PATH_MOVE_TO) {
            segment->length_to = prev_length;
        } else {
            segment->length_to = prev_length + sqrt(dx * dx + dy * dy);
        }
    }
}

ok_motion_path_t *ok_motion_path_create(const ok_path_t *path) {
    size_t first_subpath = 0;
    size_t last_subpath = path->subpaths.length - 1;
    size_t count = _ok_path_flatten_generic(path, first_subpath, last_subpath, true, false,
                                            NULL, NULL);
    ok_motion_path_t *out_path = calloc(1, sizeof(ok_motion_path_t));
    if (!vector_ensure_capacity(&out_path->segments, count)) {
        ok_motion_path_free(out_path);
        return NULL;
    } else {
        _ok_path_flatten_generic(path, first_subpath, last_subpath, true, false,
                                 _ok_motion_path_add_point, out_path);
        return out_path;
    }
}

double ok_motion_path_length(const ok_motion_path_t *path) {
    if (path->segments.length > 0) {
        return vector_last(&path->segments)->length_to;
    } else {
        return 0.0;
    }
}

static double _ok_normalize_angle(double radians) {
    if (radians < -M_PI || radians > M_PI) {
        return radians - (2.0 * M_PI) * floor((radians + M_PI) / (2.0 * M_PI));
    } else {
        return radians;
    }
}

void ok_motion_path_location(const ok_motion_path_t *path, double p,
                             double *out_x, double *out_y, double *out_angle) {
    const size_t count = path->segments.length;
    if (count == 0) {
        if (out_x) {
            *out_x = 0.0;
        }
        if (out_y) {
            *out_y = 0.0;
        }
        if (out_angle) {
            *out_angle = 0.0;
        }
    } else {
        struct ok_motion_path_segment *segments = path->segments.values;
        const double length = segments[count - 1].length_to;
        size_t p_low;
        size_t p_high;
        if (p <= 0.0) {
            p_low = 0;
            p_high = MIN(1, count - 1);
        } else if (p >= 1.0) {
            p_low = count - 1;
            p_high = count - 1;
        } else {
            // Binary search
            const double l = p * length;
            p_low = 0;
            p_high = count - 1;
            while (p_high - p_low > 1) {
                size_t index = (p_high + p_low) / 2;
                if (segments[index].length_to > l) {
                    p_high = index;
                } else {
                    p_low = index;
                }
            }
        }

        struct ok_motion_path_segment *s1 = &segments[p_low];
        struct ok_motion_path_segment *s2 = &segments[p_high];
        if (p_low == p_high || length == 0.0) {
            if (out_x) {
                *out_x = s1->x;
            }
            if (out_y) {
                *out_y = s1->y;
            }
            if (out_angle) {
                *out_angle = s2->angle_to;
            }
        } else {
            const double p1 = s1->length_to / length;
            const double p2 = s2->length_to / length;

            if (p1 == p2) {
                if (out_x) {
                    *out_x = s1->x;
                }
                if (out_y) {
                    *out_y = s1->y;
                }
                if (out_angle) {
                    *out_angle = s2->angle_to;
                }
            } else {
                if (out_x) {
                    *out_x = s1->x + (s2->x - s1->x) * (p - p1) / (p2 - p1);
                }
                if (out_y) {
                    *out_y = s1->y + (s2->y - s1->y) * (p - p1) / (p2 - p1);
                }
                if (out_angle) {
                    if (s1->type == OK_PATH_MOVE_TO || !ok_path_type_is_curve(s2->type)) {
                        *out_angle = s2->angle_to;
                    } else {
                        const double angle1 = s1->angle_to;
                        const double angle2 = s2->angle_to;
                        const double d_angle = _ok_normalize_angle(angle2 - angle1);
                        *out_angle = angle1 + d_angle * (p - p1) / (p2 - p1);
                    }
                }
            }
        }
    }
}

// MARK: Planar straight-line graph

struct ok_pslg_context {
    float *points;
    size_t num_points;
    size_t max_points;

    void *segments;
    size_t num_segments;
    size_t max_segments;

    size_t index_size;
};

static bool _ok_index_size_valid(size_t index_size) {
    // On 32-bit targets, index_size must be 1, 2, or 4.
    // On 64-bit targets, index_size must be 1, 2, 4, or 8.
    size_t test_index_size = 1;
    while (test_index_size <= sizeof(size_t)) {
        if (index_size == test_index_size) {
            return true;
        } else {
            test_index_size <<= 1;
        }
    }
    return false;
}

static void _ok_pslg_add_point(enum ok_path_element_type type, double x, double y, void *userData) {
    struct ok_pslg_context *context = userData;
    if (context->num_points < context->max_points) {
        // Add point
        size_t i = context->num_points * 2;
        context->points[i] = (float)x;
        context->points[i + 1] = (float)y;

        // Add segment
        if (type != OK_PATH_MOVE_TO && context->num_points > 0 &&
            context->num_segments < context->max_segments) {
            i = context->num_segments * 2;
            switch (context->index_size) {
                case 1: {
                    uint8_t *segments = context->segments;
                    segments[i] = (uint8_t)context->num_points - 1;
                    segments[i + 1] = (uint8_t)context->num_points;
                    break;
                }
                case 2: {
                    uint16_t *segments = context->segments;
                    segments[i] = (uint16_t)context->num_points - 1;
                    segments[i + 1] = (uint16_t)context->num_points;
                    break;
                }
                case 4: {
                    uint32_t *segments = context->segments;
                    segments[i] = (uint32_t)context->num_points - 1;
                    segments[i + 1] = (uint32_t)context->num_points;
                    break;
                }
                case 8: {
                    uint64_t *segments = context->segments;
                    segments[i] = (uint64_t)context->num_points - 1;
                    segments[i + 1] = (uint64_t)context->num_points;
                    break;
                }
            }
            context->num_segments++;
        }

        context->num_points++;
    }
}

void ok_path_create_pslg_generic(const ok_path_t *path, bool close_subpaths, size_t index_size,
                                 float **out_points, size_t *out_num_points,
                                 void **out_segment_indices, size_t *out_num_segments) {
    // Validate input
    if (!_ok_index_size_valid(index_size) || path->elements.length == 0) {
        *out_points = NULL;
        *out_segment_indices = NULL;
        *out_num_points = 0;
        *out_num_segments = 0;
        return;
    }

    size_t first_subpath = 0;
    size_t last_subpath = path->subpaths.length - 1;

    // Initialize context
    struct ok_pslg_context context;
    context.index_size = index_size;
    context.points = NULL;
    context.segments = NULL;
    context.num_points = 0;
    context.num_segments = 0;
    context.max_points = _ok_path_flatten_generic(path, first_subpath, last_subpath,
                                                  true, close_subpaths, NULL, NULL);
    context.max_segments = MIN((context.max_points - 1), (((size_t)1) << (index_size * 8 - 1)));

    // Allocate
    context.points = malloc(context.max_points * sizeof(float[2]));
    context.segments = malloc(context.max_segments * 2 * index_size);
    if (!context.points || !context.segments) {
        free(context.points);
        free(context.segments);
        *out_points = NULL;
        *out_segment_indices = NULL;
        *out_num_points = 0;
        *out_num_segments = 0;
        return;
    }

    // Convert
    _ok_path_flatten_generic(path, first_subpath, last_subpath, true, close_subpaths,
                             _ok_pslg_add_point, &context);
    *out_points = context.points;
    *out_num_points = context.num_points;
    *out_segment_indices = context.segments;
    *out_num_segments = context.num_segments;
}

void ok_path_append_pslg_generic(ok_path_t *path, size_t index_size, const float *points,
                                 const void *segments, size_t num_segments) {
    if (!_ok_index_size_valid(index_size) || num_segments == 0) {
        return;
    }
    size_t last_segment = (size_t)-1;
    for (size_t i = 0; i < num_segments; i++) {
        size_t p1 = 0;
        size_t p2 = 0;
        switch (index_size) {
            case 1: {
                const uint8_t *segments8 = segments;
                p1 = segments8[i * 2];
                p2 = segments8[i * 2 + 1];
                break;
            }
            case 2: {
                const uint16_t *segments16 = segments;
                p1 = segments16[i * 2];
                p2 = segments16[i * 2 + 1];
                break;
            }
            case 4: {
                const uint32_t *segments32 = segments;
                p1 = segments32[i * 2];
                p2 = segments32[i * 2 + 1];
                break;
            }
            case 8: {
                const uint64_t *segments64 = segments;
                p1 = (size_t)segments64[i * 2];
                p2 = (size_t)segments64[i * 2 + 1];
                break;
            }
        }

        if (p1 != last_segment) {
            ok_path_move_to(path, (double)points[p1 * 2], (double)points[p1 * 2 + 1]);
        }
        ok_path_line_to(path, (double)points[p2 * 2], (double)points[p2 * 2 + 1]);

        last_segment = p2;
    }
}
