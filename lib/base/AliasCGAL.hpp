/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

// This file is to collect all CGAL aliases inside yade sub-namespace

#ifdef YADE_CGAL

#include <lib/high-precision/CgalNumTraits.hpp>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/squared_distance_3.h>

namespace yade {

// CGAL definitions - does not work with another kernel!! Why???
//   / Answer: because to use exact kernal yade must link with mpfr / Janek.
//   / to do this have a look at CMakeLists 'Real precision' section

	// These are taken from files pkg/dem/Polyhedra.hpp and /pkg/dem/Gl1_PotentialBlock.hpp
	// let's call this CGP for CGal Polyhedrons namespace. It uses inexact kernel. Later when mpfr support is added we can have similar
	// aliases in other namespace for exact kernel. Or maybe later rename namespace into struct, so that this can be used as template argument.
	// We will see how this will need to evelove.
//	namespace CGP {
		using K              = CGAL::Exact_Real_predicates_inexact_constructions_kernel;
		using Polyhedron     = CGAL::Polyhedron_3<K>;
		using Mesh           = CGAL::Surface_mesh<K::Point_3>;
		using Triangulation  = CGAL::Delaunay_triangulation_3<K>;
		using CGALpoint      = K::Point_3;
		using CGALtriangle   = K::Triangle_3;
		using CGALvector     = K::Vector_3;
		using Transformation = CGAL::Aff_transformation_3<K>;
		using Segment        = K::Segment_3;
		using Triangle       = CGAL::Triangle_3<K>;
		using Plane          = CGAL::Plane_3<K>;
		using Line           = CGAL::Line_3<K>;
		using CGAL_ORIGIN    = CGAL::Origin;
		using CGAL_AABB_tree = CGAL::AABB_tree<CGAL::AABB_traits<K,CGAL::AABB_triangle_primitive<K,std::vector<Triangle>::iterator>>>;
//	}
	// It would be perhaps useful to collect here similar stuff frmo other files.
} // namespace yade

#endif

