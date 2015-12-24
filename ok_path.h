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

#ifndef _OK_PATH_H_
#define _OK_PATH_H_

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif
    
    typedef struct ok_path ok_path;
    
    ok_path *ok_path_alloc(void);
    
    void ok_path_free(ok_path *path);
    
    bool ok_path_equals(const ok_path *path1, const ok_path *path2);
    
    // MARK: Path creation
    
    void ok_path_move_to(ok_path *path, const double x, const double y);
    
    void ok_path_line_to(ok_path *path, const double x, const double y);
    
    void ok_path_curve_to(ok_path *path, const double cx1, const double cy1, const double cx2, const double cy2,
                          const double x, const double y);
    
    void ok_path_quad_curve_to(ok_path *path, const double cx, const double cy, const double x, const double y);
    
    void ok_path_elliptical_arc_to(ok_path *path, const double radius_x, const double radius_y,
                                   const double rotation_radians, const bool large_arc, const bool sweep,
                                   const double x, const double y);
    
    void ok_path_append(ok_path *path, ok_path *path_to_append);
    
    /** 
     Append an SVG path, as defined in http://www.w3.org/TR/SVG/paths.html
     If successful, returns true.
     If unsuccessful, returns false and [error], if not NULL, is set to an error string.
     The error string is maintained internally and should not be freed.
     */
    bool ok_path_append_svg(ok_path *path, const char *svg_path, char **error);

    void ok_path_close(ok_path *path);
    
    // MARK: Path querying and flattening
    
    /// Gets the length of the path. First, the path is flattened if it hasn't been flattened yet.
    double ok_path_get_length(ok_path *path);
    
#ifdef __cplusplus
}
#endif

#endif
