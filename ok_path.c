/*
 ok-path
 Copyright (c) 2016 David Brackeen
 
 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from the
 use of this software. Permission is granted to anyone to use this software
 for any purpose, including commercial applications, and to alter it and
 redistribute it freely, subject to the following restrictions:
 
 1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software in a
    product, an acknowledgment in the product documentation would be appreciated
    but is not required.
 2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.
 3. This notice may not be removed or altered from any source distribution.
 */

#include "ok_path.h"
#include <ctype.h>
#include <math.h>
#include <memory.h>
#include <stdio.h> // For snprintf
#include <stdlib.h>

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

// MARK: vector

typedef struct {
    void *values;
    unsigned int length;
    unsigned int capacity;
} ok_vector;

static bool ok_vector_ensure_capacity(ok_vector *list, const size_t element_size, const unsigned int additional_count) {
    if (list) {
        if (!list->values || list->length + additional_count > list->capacity) {
            const unsigned int new_capacity = MAX(16, MAX(list->length + additional_count, list->capacity << 1));
            void *new_data = realloc(list->values, element_size * new_capacity);
            if (!new_data) {
                return false;
            }
            list->values = new_data;
            list->capacity = new_capacity;
        }
        return true;
    }
    return false;
}

static void ok_vector_free(ok_vector *list) {
    if (list && list->values) {
        free(list->values);
        list->values = NULL;
        list->length = 0;
        list->capacity = 0;
    }
}

// MARK: path

typedef enum {
    MOVE_TO = 0,
    LINE_TO,
    CURVE_TO,
} ok_path_type;

typedef struct {
    ok_path_type type;
    double x, y;
    
    // Control points. Valid if type is CURVE_TO, otherwise these are undefined.
    double cx1, cy1;
    double cx2, cy2;
} ok_path_segment;

struct ok_path {
    ok_vector path_segments;
    double subpath_origin_x;
    double subpath_origin_y;
};

ok_path *ok_path_alloc() {
    return calloc(1, sizeof(ok_path));
}

void ok_path_free(ok_path *path) {
    ok_vector_free(&path->path_segments);
    free(path);
}

bool ok_path_equals(const ok_path *path1, const ok_path *path2) {
    if (path1->path_segments.length != path2->path_segments.length) {
        return false;
    }
    ok_path_segment *segment1 = path1->path_segments.values;
    ok_path_segment *segment2 = path2->path_segments.values;
    for (unsigned int i = 0; i < path1->path_segments.length; i++) {
        if (segment1->type != segment2->type) {
            return false;
        }
        if (segment1->x != segment2->x || segment1->y != segment2->y) {
            return false;
        }
        if (segment1->type == CURVE_TO) {
            if (segment1->cx1 != segment2->cx1 || segment1->cy1 != segment2->cy1 ||
                segment1->cx2 != segment2->cx2 || segment1->cy2 != segment2->cy2) {
                return false;
            }
        }
        segment1++;
        segment2++;
    }
    return true;
}

static double ok_path_last_x(const ok_path *path) {
    if (path->path_segments.length) {
        ok_path_segment *values = path->path_segments.values;
        return values[path->path_segments.length - 1].x;
    } else {
        return 0;
    }
}

static double ok_path_last_y(const ok_path *path) {
    if (path->path_segments.length) {
        ok_path_segment *values = path->path_segments.values;
        return values[path->path_segments.length - 1].y;
    } else {
        return 0;
    }
}

// MARK: Commands

void ok_path_move_to(ok_path *path, const double x, const double y) {
    if (ok_vector_ensure_capacity(&path->path_segments, sizeof(ok_path_segment), 1)) {
        ok_path_segment *values = path->path_segments.values;
        ok_path_segment *segment = &values[path->path_segments.length++];
        segment->type = MOVE_TO;
        segment->x = x;
        segment->y = y;
        path->subpath_origin_x = x;
        path->subpath_origin_y = y;
    }
}

void ok_path_line_to(ok_path *path, const double x, const double y) {
    if (ok_vector_ensure_capacity(&path->path_segments, sizeof(ok_path_segment), 1)) {
        ok_path_segment *values = path->path_segments.values;
        ok_path_segment *segment = &values[path->path_segments.length++];
        segment->type = LINE_TO;
        segment->x = x;
        segment->y = y;
    }
}

void ok_path_curve_to(ok_path *path, const double cx1, const double cy1, const double cx2, const double cy2,
                      const double x, const double y) {
    if (ok_vector_ensure_capacity(&path->path_segments, sizeof(ok_path_segment), 1)) {
        ok_path_segment *values = path->path_segments.values;
        ok_path_segment *segment = &values[path->path_segments.length++];
        segment->type = CURVE_TO;
        segment->x = x;
        segment->y = y;
        segment->cx1 = cx1;
        segment->cy1 = cy1;
        segment->cx2 = cx2;
        segment->cy2 = cy2;
    }
}

void ok_path_append(ok_path *path, ok_path *path_to_append) {
    unsigned int count = path_to_append->path_segments.length;
    if (ok_vector_ensure_capacity(&path->path_segments, sizeof(ok_path_segment), count)) {
        ok_path_segment *values = path->path_segments.values;
        ok_path_segment *src_values = path_to_append->path_segments.values;
        ok_path_segment *dst_values = &values[path->path_segments.length];
        memcpy(dst_values, src_values, count * sizeof(ok_path_segment));
        path->path_segments.length += count;
        path->subpath_origin_x = path_to_append->subpath_origin_x;
        path->subpath_origin_y = path_to_append->subpath_origin_y;
    }
}

void ok_path_close(ok_path *path) {
    ok_path_line_to(path, path->subpath_origin_x, path->subpath_origin_y);
}

void ok_path_quad_curve_to(ok_path *path, const double cx, const double cy, const double x, const double y) {
    double last_x = ok_path_last_x(path);
    double last_y = ok_path_last_y(path);
    double cx1 = (last_x + cx * 2.0) / 3.0;
    double cy1 = (last_y + cy * 2.0) / 3.0;
    double cx2 = (x + cx * 2.0) / 3.0;
    double cy2 = (y + cy * 2.0) / 3.0;
    ok_path_curve_to(path, cx1, cy1, cx2, cy2, x, y);
}

void ok_path_elliptical_arc_to(ok_path *path, const double radius_x, const double radius_y,
                               const double rotation_radians, const bool large_arc, const bool sweep,
                               const double x, const double y) {
    // Convert arc to a series of bezier curves.
    // Interface is the same as the SVG spec, except xAxisRotation is in radians, not degrees.
    // See http://www.w3.org/TR/SVG/paths.html#PathDataEllipticalArcCommands
    // See http://www.w3.org/TR/SVG/implnote.html#ArcImplementationNotes
    if (radius_x == 0 || radius_y == 0) {
        ok_path_line_to(path, x, y);
    } else {
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
        if (uy < 0) {
            angle_start = -angle_start;
        }
        n = sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy));
        p = ux * vx + uy * vy;
        double angle_extent = acos(p / n);
        if (ux * vy - uy * vx < 0) {
            angle_extent = -angle_extent;
        }
        if (!sweep && angle_extent > 0) {
            angle_extent -= 2.0 * M_PI;
        } else if (sweep && angle_extent < 0) {
            angle_extent += 2.0 * M_PI;
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
        
        double s = angle_start;
        double d = (sweep ? M_PI_2 : -M_PI_2);
        while (true) {
            if ((sweep && s > angle_end) || (!sweep && s < angle_end)) {
                break;
            }
            
            double e = s + d;
            if ((sweep && e > angle_end) || (!sweep && e < angle_end)) {
                e = angle_end;
            }
            
            double da = e - s;
            double alpha = 4.0 * tan(da / 4.0) / 3.0;
            
            // Alternative alpha from: http://www.spaceroots.org/documents/ellipse/
            // TODO: test for very large arcs to see which alpha is better
            //double t = tan(da / 2.0);
            //double alpha = sin(da) * (fsqrt(4.0 + 3.0 * t * t) - 1.0) / 3.0;
            
            double x_a = x_b;
            double y_a = y_b;
            double x_a_dot = x_b_dot;
            double y_a_dot = y_b_dot;
            
            cos_eta_b = cos(e);
            sin_eta_b = sin(e);
            
            a_cos_eta_b = rx * cos_eta_b;
            b_sin_eta_b = ry * sin_eta_b;
            a_sin_eta_b = rx * sin_eta_b;
            b_cos_eta_b = ry * cos_eta_b;
            
            x_b_dot = -a_sin_eta_b * cos_a - b_cos_eta_b * sin_a;
            y_b_dot = -a_sin_eta_b * sin_a + b_cos_eta_b * cos_a;
            
            if (e == angle_end) {
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
            
            s += d;
        }
    }
}

// MARK: SVG parsing

static const char * const SVG_COMMANDS = "MmLlHhVvAaQqTtCcSsZz";
static const unsigned int SVG_COMMAND_VALUES[] = { 2, 2, 2, 2, 1, 1, 1, 1, 7, 7, 4, 4, 2, 2, 6, 6, 4, 4, 0, 0 };
static const char * const SVG_WHITESPACE = " \t\r\n";
static const char * const SVG_WHITESPACE_OR_COMMA = ", \t\r\n";

static char SVG_ERROR_MESSAGE[80];

bool ok_path_append_svg(ok_path *path, const char *svg_path, char **error) {
    const char *str = svg_path;
    double values[7];
    char last_command = 0;
    unsigned int last_values_required = 0;
    bool last_control_set = false;
    double curr_x = 0, curr_y = 0;
    double control_x = 0, control_y = 0;
    
    if (!str && error) {
        *error = "NULL input";
    }
    
    while (str && *str) {
        bool command_required = last_values_required == 0;
        const char *skip_chars = command_required ? SVG_WHITESPACE : SVG_WHITESPACE_OR_COMMA;
        if (strchr(skip_chars, *str)) {
            str++;
        } else {
            // Get command
            char command;
            unsigned int values_required;
            char *found = strchr(SVG_COMMANDS, *str);
            if (found || command_required) {
                command = *str++;
                values_required = found ? SVG_COMMAND_VALUES[found - SVG_COMMANDS] : 0;
            } else {
                command = last_command;
                values_required = last_values_required;
            }
            
            // Parse numbers
            if (values_required > 0) {
                unsigned int count = 0;
                while (*str && count < values_required) {
                    if (count > 0 && strchr(SVG_WHITESPACE_OR_COMMA, *str)) {
                        str++;
                    } else {
                        char *endptr;
                        values[count++] = strtod(str, &endptr);
                        if (str == endptr) {
                            if (error) {
                                snprintf(SVG_ERROR_MESSAGE, sizeof(SVG_ERROR_MESSAGE),
                                         "Could not parse number at position: %li", str - svg_path);
                                *error = SVG_ERROR_MESSAGE;
                            }
                            return false;
                        }
                        str = endptr;
                    }
                }
                if (count < values_required) {
                    if (error) {
                        snprintf(SVG_ERROR_MESSAGE, sizeof(SVG_ERROR_MESSAGE),
                                 "Unexpected EOF: Needed %i more numbers", (values_required - count));
                        *error = SVG_ERROR_MESSAGE;
                    }
                    return false;
                }
            }
            
            // Execute command
            bool control_set = false;
            switch (command) {
                case 'Z': case 'z': {
                    curr_x = path->subpath_origin_x;
                    curr_y = path->subpath_origin_y;
                    ok_path_line_to(path, curr_x, curr_y);
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
                case 'S': case 's': {
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
                case 'T': case 't': {
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
                case 'A': case 'a': {
                    double radius_x = values[0];
                    double radius_y = values[1];
                    double angle = values[2] * M_PI / 180.0;
                    bool large_arc = values[3] != 0;
                    bool sweep = values[4] != 0;
                    if (command == 'A') {
                        curr_x = values[5];
                        curr_y = values[6];
                    } else {
                        curr_x += values[5];
                        curr_y += values[6];
                    }
                    ok_path_elliptical_arc_to(path, radius_x, radius_y, angle, large_arc, sweep, curr_x, curr_y);
                    break;
                }
                default: {
                    if (error) {
                        snprintf(SVG_ERROR_MESSAGE, sizeof(SVG_ERROR_MESSAGE),
                                 "Invalid SVG command %c at position %li", command, (str-svg_path));
                        *error = SVG_ERROR_MESSAGE;
                    }
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
