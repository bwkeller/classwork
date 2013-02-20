/*============================================================================
 * Define (conforming or non-conforming) mesh joinings.
 *============================================================================*/

/* Code_Saturne version 2.1.4 */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_join.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define mesh joinings.
 *
 * This is done by calling the cs_join_add() function for each joining
 * operation to add.
 *
 * The arguments to cs_join_add() are:
 *   sel_criteria <-- boundary face selection criteria string
 *   fraction     <-- value of the fraction parameter;
 *                    the initial tolerance radius associated to each vertex
 *                    is equal to the lenght of the shortest incident edge,
 *                    multiplied by this fraction.
 *   plane        <-- value of the plane parameter;
 *                    when subdividing faces, 2 faces are considered
 *                    coplanar and may be joined if angle between their
 *                    normals (in degrees) does not exceed this parameter.
 *   verbosity    <-- level of verbosity required
 *
 * The function returns a number (1 to n) associated with the
 * new joining. This number may optionnally be used to assign advanced
 * parameters to the joining.
 *----------------------------------------------------------------------------*/

void
cs_user_join(void)
{
  int    join_num;
  int    verbosity = 1, visualization;
  float  fraction = 0.10, plane = 0.25;

  fraction = 0.10;
  plane = 25.0;
  verbosity = 1; /* debug level if >= 3 */
  visualization = 1; /* debug level if >= 3 */

  /* Add a joining operation */
  /* ----------------------- */

  join_num = cs_join_add("98 or 99",
                         fraction,
                         plane,
                         verbosity,
                         visualization);

  /*--------------------------------------------------------------------------*/

  /* Example with advanced parameters;
     Advanced parameters may be modified to solve errors during the
     joining step or to get a better mesh quality. */

  {
    /* Merge tolerance factor:
     * used to locally modify the tolerance associated to each
     * vertex BEFORE the merge step.
     *   = 0 => no vertex merge;
     *   < 1 => vertex merge is more strict. It may increase the number
     *          of tolerance reduction and so define smaller subset of
     *          vertices to merge together but it can drive to loose
     *          intersections;
     *   = 1 => no change;
     *   > 1 => vertex merge is less strict. The subset of vertices able
     *          to be merged together is greater. */

     double mtf = 1.00;

     /* Pre-merge factor:
      * this parameter is used to define a bound under which two vertices
      * are merged before the merge step.
      * Tolerance limit for the pre-merge = pmf * fraction. */

     double pmf = 0.10;

     /* Tolerance computation mode:
      *
      *   1: (default) tol = min. edge length related to a vertex * fraction
      *   2: tolerance is computed like in mode 1 with in addition, the
      *      multiplication by a coefficient equal to the max sin(e1, e2)
      *      where e1 and e2 are two edges sharing the same vertex V for which
      *      we want to compute the tolerance.
      *  11: as 1 but taking into account only the selected faces
      *  12: as 2 but taking into account only the selected faces */

     int tcm = 1;

      /* Intersection computation mode:
       *  1: (default) Original algorithm.
       *     Try to snap intersection to extremity.
       *  2: New intersection algorithm.
       *     Avoid snapping intersection to extremity. */

     int icm = 1;

     /* Maximum number of equivalence breaks which is
      * enabled during the merge step. */

     int max_break = 500;

     /* Maximum number of sub-faces when splitting a selected face. */

     int max_sub_face = 100;

     /* tml, tmb and tmr are parameters of the searching algorithm for
      * face intersections between selected faces (octree structure).
      * Useful to adjust speed vs. memory consumption. */

     /* Tree Max Level:
      * deepest level reachable when building the tree */

     int tml = 30;

     /* Tree Max Boxes:
      * max. number of bounding boxes which can be linked to a leaf of the tree
      * (not necessary true for the deepest level) */

     int tmb = 25;

     /* Tree Max. Ratio:
      * stop refining the tree structure when
      * number of bounding boxes > tmr * number of faces to locate.
      * In parallel, a separate (usually lower) value may be set for
      * the initial coarse tree used to determine distribution.
      * Reducing this will help reduce memory consumption. */

     double tmr = 5.0;
     double tmr_distrib = 2.0;

     /* Set advanced parameters */

     cs_join_set_advanced_param(join_num,
                                mtf, pmf, tcm, icm,
                                max_break, max_sub_face,
                                tml, tmb, tmr, tmr_distrib);
  }
}
