/*
 ok-path
 https://github.com/brackeen/ok-path
 Copyright (c) 2016-2018 David Brackeen

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

#ifndef OK_PATH_H
#define OK_PATH_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum ok_path_element_type {
    OK_PATH_MOVE_TO = 0,
    OK_PATH_LINE_TO,
    OK_PATH_QUAD_CURVE_TO,
    OK_PATH_CUBIC_CURVE_TO,
    OK_PATH_CLOSE
};

/**
 * An `ok_path_t` is a series of line segments and Bézier curves.
 */
typedef struct ok_path ok_path_t;

typedef struct ok_motion_path ok_motion_path_t;

// MARK: Creating paths

/**
 * Creates a new #ok_path_t.
 * The path should be freed with #ok_path_free(ok_path_t *).
 */
ok_path_t *ok_path_create(void);

/**
 * Frees an #ok_path_t.
 */
void ok_path_free(ok_path_t *path);

// MARK: Modifying paths

/**
 * Moves the current location of the path, creating a new subpath.
 *
 * @param path The path.
 * @param x The x location to move to.
 * @param y The y location to move to.
 */
void ok_path_move_to(ok_path_t *path, double x, double y);

/**
 * Adds a line segment to the path.
 *
 * @param path The path.
 * @param x The x location of the end point of the line segment.
 * @param y The y location of the end point of the line segment.
 */
void ok_path_line_to(ok_path_t *path, double x, double y);

/**
 * Adds a cubic Bézier curve to the path.
 *
 * @param path The path.
 * @param cx1 The x location of the first control point.
 * @param cy1 The y location of the first control point.
 * @param cx2 The x location of the second control point.
 * @param cy2 The y location of the second control point.
 * @param x The x location of the end point of the curve.
 * @param y The y location of the end point of the curve.
 */
void ok_path_curve_to(ok_path_t *path, double cx1, double cy1, double cx2, double cy2,
                      double x, double y);

/**
 * Adds a quadratic Bézier curve to the path.
 *
 * @param path The path.
 * @param cx The x location of the control point.
 * @param cy The y location of the control point.
 * @param x The x location of the end point of the curve.
 * @param y The y location of the end point of the curve.
 */
void ok_path_quad_curve_to(ok_path_t *path, double cx, double cy, double x, double y);

/**
 * Adds an arc to the path.
 *
 * @param path The path.
 * @param radius The radius of the arc.
 * @param large_arc A flag indicating that the larger arc (greater than 180 degrees) is drawn.
 * @param sweep A flag indicating that the arc is drawn in the "positive angle" direction.
 * @param x The x location of the end point of the arc.
 * @param y The y location of the end point of the arc.
 * @see http://www.w3.org/TR/SVG/paths.html#PathDataEllipticalArcCommands
 */
void ok_path_arc_to(ok_path_t *path, double radius, bool large_arc, bool sweep, double x, double y);

/**
 * Adds an elliptical arc to the path.
 *
 * @param path The path.
 * @param radius_x The radius on the x axis.
 * @param radius_y The radius on the y axis.
 * @param rotation_radians The angle, in radians, of the arc.
 * @param large_arc A flag indicating that the larger arc (greater than 180 degrees) is drawn.
 * @param sweep A flag indicating that the arc is drawn in the "positive angle" direction.
 * @param x The x location of the end point of the arc.
 * @param y The y location of the end point of the arc.
 * @see http://www.w3.org/TR/SVG/paths.html#PathDataEllipticalArcCommands
 */
void ok_path_elliptical_arc_to(ok_path_t *path, double radius_x, double radius_y,
                               double rotation_radians, bool large_arc, bool sweep,
                               double x, double y);

/**
 * Appends a path to the path.
 * @param path The path.
 * @param path_to_append The path to append.
 */
void ok_path_append(ok_path_t *path, const ok_path_t *path_to_append);

/**
 * Appends a series of line segments to the path. Example:
 *
 *     double points[][2] = {{0, 0}, {10, 20}, {40, 20}, {30, 10}, {50, 0}, {0, 0}};
 *     ok_path_append_lines(path, OK_PATH_MOVE_TO, points, sizeof(points) / sizeof(*points));
 *
 * @param path The path.
 * @param first_type The element type of the first point, which can be either #OK_PATH_MOVE_TO
 *     or #OK_PATH_LINE_TO.
 * @param points The array of points.
 * @param num_points The number of points.
 */
void ok_path_append_lines(ok_path_t *path, enum ok_path_element_type first_type,
                          const double (*points)[2], size_t num_points);

/**
 * Appends an SVG path, like "M 100,100 L 200,100 200,200 100,200 Z".
 * The entire specification defined at http://www.w3.org/TR/SVG/paths.html is accepted.
 *
 * If successful, returns `true`. Otherwise, returns `false` and `out_error_message` is set to an
 * error string. The error string is maintained internally and should not be freed.
 *
 * @param path The path.
 * @param svg_path The SVG path string to parse.
 * @param[out] out_error_message The pointer where the error string should be stored. May be `NULL`.
 * @return `true` if successful.
 */
bool ok_path_append_svg(ok_path_t *path, const char *svg_path, char **out_error_message);

/**
 * Closes the current subpath.
 */
void ok_path_close(ok_path_t *path);

/**
 * Remove all elements from the path.
 */
void ok_path_reset(ok_path_t *path);

/**
 * Creates a flattened version of this path. All curved segments of the path are converted to a
 * series of straight lines that approximates the curve.
 */
ok_path_t *ok_path_flatten(const ok_path_t *path);

// MARK: Getting information about paths

/**
 * Checks if two paths are equal.
 *
 * @return `true` if both paths contain the same sequence of path elements.
 */
bool ok_path_equals(const ok_path_t *path1, const ok_path_t *path2);

/**
 * Checks if a path is flat (consists of only #OK_PATH_MOVE_TO, #OK_PATH_LINE_TO, and #OK_PATH_CLOSE
 * elements).
 */
bool ok_path_is_flat(const ok_path_t *path);

/**
 * Gets the number of elements in the path.
 */
size_t ok_path_element_count(const ok_path_t *path);

/**
 * Gets an element in the path.
 *
 * @param path The path.
 * @param index The element index, from `0` to `count-1`, where `count` is the value returned from
 *     #ok_path_element_count.
 * @param[out] out_cx1 The x location of the first control point. This value is only set if the
 *     element type is #OK_PATH_QUAD_CURVE_TO or #OK_PATH_CUBIC_CURVE_TO.
 * @param[out] out_cy1 The y location of the first control point. This value is only set if the
 *     element type is #OK_PATH_QUAD_CURVE_TO or #OK_PATH_CUBIC_CURVE_TO.
 * @param[out] out_cx2 The x location of the second control point. This value is only set if the
 *     element type is #OK_PATH_CUBIC_CURVE_TO.
 * @param[out] out_cy2 The y location of the second control point. This value is only set if the
 *     element type is #OK_PATH_CUBIC_CURVE_TO.
 * @param[out] out_x The x location of the element's end point.
 * @param[out] out_y The y location of the element's end point.
 * @return the element type.
 */
enum ok_path_element_type ok_path_element_get(const ok_path_t *path, size_t index,
                                              double *out_cx1, double *out_cy1,
                                              double *out_cx2, double *out_cy2,
                                              double *out_x, double *out_y);

/**
 * Sets the points associated with an element in the path.
 * @param index The element index, from `0` to `count-1`, where `count` is the value returned from
 *     #ok_path_element_count.
 * @param cx1 The x location of the first control point. This value is only used if the
 *     element type is #OK_PATH_QUAD_CURVE_TO or #OK_PATH_CUBIC_CURVE_TO.
 * @param cy1 The y location of the first control point. This value is only used if the
 *     element type is #OK_PATH_QUAD_CURVE_TO or #OK_PATH_CUBIC_CURVE_TO.
 * @param cx2 The x location of the second control point. This value is only used if the
 *     element type is #OK_PATH_CUBIC_CURVE_TO.
 * @param cy2 The y location of the second control point. This value is only used if the
 *     element type is #OK_PATH_CUBIC_CURVE_TO.
 * @param x The x location of the element's end point.
 * @param y The y location of the element's end point.
 */
void ok_path_element_set(const ok_path_t *path, size_t index,
                         double cx1, double cy1, double cx2, double cy2, double x, double y);

// MARK: Subpaths

/**
 * Gets the number of subpaths in a path.
 */
size_t ok_subpath_count(const ok_path_t *path);

/**
 * Creates a copy of a subpath from the path.
 *
 * @param path The path.
 * @param index The subpath index, from `0` to `count-1`, where `count` is the value returned from
 *     #ok_subpath_count.
 */
ok_path_t *ok_subpath_create(const ok_path_t *path, size_t index);

/**
 * Creates a flattened version of a subpath. All curved segments of the subpath are converted to a
 * series of straight lines that approximates the curve.
 *
 * @param path The path.
 * @param index The subpath index, from `0` to `count-1`, where `count` is the value returned from
 *     #ok_subpath_count.
 */
ok_path_t *ok_subpath_flatten(const ok_path_t *path, size_t index);

/**
 * Gets the index of the first element of the specified subpath in the path.
 * @param path The path.
 * @param index The subpath index, from `0` to `count-1`, where `count` is the value returned from
 *     #ok_subpath_count.
 */
size_t ok_subpath_first_element_index(const ok_path_t *path, size_t index);

/**
 * Gets the index of the last element of the specified subpath in the path.
 * @param path The path.
 * @param index The subpath index, from `0` to `count-1`, where `count` is the value returned from
 *     #ok_subpath_count.
 */
size_t ok_subpath_last_element_index(const ok_path_t *path, size_t index);

/**
 * Checks if a subpath is flat (consists of only #OK_PATH_MOVE_TO, #OK_PATH_LINE_TO, and
 * #OK_PATH_CLOSE elements).
 */
bool ok_subpath_is_flat(const ok_path_t *path, size_t index);

// MARK: Motion paths

/**
 * Creates a motion path from an existing path. A motion path can be used to get a location
 * along the path, get the angle of the path at that location, or to animate a point
 * along the path at a constant rate.
 */
ok_motion_path_t *ok_motion_path_create(const ok_path_t *path);

/**
 * Frees the path.
 */
void ok_motion_path_free(ok_motion_path_t *path);

/**
 * Gets the length of the path.
 */
double ok_motion_path_length(const ok_motion_path_t *path);

/**
 * Gets the location on the path at location `p`, where `p` is a value from 0.0 to 1.0.
 *
 * @param path The path.
 * @param p The location on the path, from 0.0 to 1.0.
 * @param[out] out_x The pointer where the x location should be stored. May be `NULL`.
 * @param[out] out_y The pointer where the y location should be stored. May be `NULL`.
 * @param[out] out_angle The pointer where the angle, in radians, should be stored. May be `NULL`.
 */
void ok_motion_path_location(const ok_motion_path_t *path, double p,
                             double *out_x, double *out_y, double *out_angle);

#ifdef __cplusplus
}
#endif

#endif
