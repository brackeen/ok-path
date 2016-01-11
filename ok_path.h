/*
 ok-path
 https://github.com/brackeen/ok-path
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

#ifndef _OK_PATH_H_
#define _OK_PATH_H_

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * An `ok_path` is a series of line segments and Bézier curves.
 */
typedef struct ok_path ok_path;

// MARK: Creating paths

/**
 * Creates a new #ok_path.
 * The path should be freed with #ok_path_free(ok_path *).
 */
ok_path *ok_path_alloc(void);

/**
 * Frees an #ok_path.
 */
void ok_path_free(ok_path *path);

// MARK: Modifying paths

/**
 * Moves the current location of the path, creating a new subpath.
 *
 * @param path The path.
 * @param x The x location to move to.
 * @param y The y location to move to.
 */
void ok_path_move_to(ok_path *path, double x, double y);

/**
 * Adds a line segment to the path.
 *
 * @param path The path.
 * @param x The x location of the end point of the line segment.
 * @param y The y location of the end point of the line segment.
 */
void ok_path_line_to(ok_path *path, double x, double y);

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
void ok_path_curve_to(ok_path *path, double cx1, double cy1, double cx2, double cy2,
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
void ok_path_quad_curve_to(ok_path *path, double cx, double cy, double x, double y);

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
void ok_path_arc_to(ok_path *path, double radius, bool large_arc, bool sweep, double x, double y);

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
void ok_path_elliptical_arc_to(ok_path *path, double radius_x, double radius_y,
                               double rotation_radians, bool large_arc, bool sweep,
                               double x, double y);

/**
 * Appends a path to the path.
 * @param path The path.
 * @param path_to_append The path to append.
 */
void ok_path_append(ok_path *path, const ok_path *path_to_append);

/**
 * Appends an SVG path, like "M 100,100 L 200,100 200,200 100,200 Z".
 * The entire specification defined at http://www.w3.org/TR/SVG/paths.html is accepted.
 *
 * If successful, returns `true`. Otherwise, returns `false` and `out_error_message` is set to an
 * error string.
 * The error string is maintained internally and should not be freed.
 *
 * @param path The path.
 * @param svg_path The SVG path string to parse.
 * @param[out] out_error_message The pointer where the error string should be stored. May be `NULL`.
 * @return `true` if successful.
 */
bool ok_path_append_svg(ok_path *path, const char *svg_path, char **out_error_message);

/**
 * Closes the current subpath.
 */
void ok_path_close(ok_path *path);

// MARK: Getting information about paths

/**
 * Checks if two paths are equal.
 *
 * @return `true` if both paths contain the same sequence of path segments.
 */
bool ok_path_equals(const ok_path *path1, const ok_path *path2);

/**
 * Gets the length of the path.
 * The path is flattened if it has not been flattened yet.
 */
double ok_path_get_length(ok_path *path);

/**
 * Gets the location on the path at location `p`, where `p` is a value from 0.0 to 1.0.
 * The path is flattened if it has not been flattened yet.
 *
 * @param path The path.
 * @param p The location on the path, from 0.0 to 1.0.
 * @param[out] out_x The pointer where the x location should be stored. May be `NULL`.
 * @param[out] out_y The pointer where the y location should be stored. May be `NULL`.
 * @param[out] out_angle The pointer where the angle, in radians, should be stored. May be `NULL`.
 */
void ok_path_get_location(ok_path *path, double p, double *out_x, double *out_y, double *out_angle);

#ifdef __cplusplus
}
#endif

#endif
