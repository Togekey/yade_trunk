/*************************************************************************
*  Copyright (C) 2008 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <lib/serialization/Serializable.hpp>
#include <pkg/common/Dispatching.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/common/Wall.hpp>

namespace yade { // Cannot have #include directive inside.

class Ig2_Facet_Sphere_ScGeom : public IGeomFunctor {
public:
	virtual bool
	             go(const shared_ptr<Shape>&       cm1,
	                const shared_ptr<Shape>&       cm2,
	                const State&                   state1,
	                const State&                   state2,
	                const Vector3r&                shift2,
	                const bool&                    force,
	                const shared_ptr<Interaction>& c);
	virtual bool goReverse(
	        const shared_ptr<Shape>&       cm1,
	        const shared_ptr<Shape>&       cm2,
	        const State&                   state1,
	        const State&                   state2,
	        const Vector3r&                shift2,
	        const bool&                    force,
	        const shared_ptr<Interaction>& c);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Ig2_Facet_Sphere_ScGeom,IGeomFunctor,"Create/update a :yref:`ScGeom` instance representing intersection of :yref:`Facet` and :yref:`Sphere`.",
		((Real,shrinkFactor,((void)"no shrinking",0),,"The radius of the inscribed circle of the facet is decreased by the value of the sphere's radius multipled by *shrinkFactor*. From the definition of contact point on the surface made of facets, the given surface is not continuous and becomes in effect surface covered with triangular tiles, with gap between the separate tiles equal to the sphere's radius multiplied by 2×*shrinkFactor*. If zero, no shrinking is done."))
	);
	// clang-format on
	DECLARE_LOGGER;
	FUNCTOR2D(Facet, Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Facet, Sphere);
};

REGISTER_SERIALIZABLE(Ig2_Facet_Sphere_ScGeom);

class Ig2_Facet_Sphere_ScGeom6D : public Ig2_Facet_Sphere_ScGeom {
public:
	virtual bool
	             go(const shared_ptr<Shape>&       cm1,
	                const shared_ptr<Shape>&       cm2,
	                const State&                   state1,
	                const State&                   state2,
	                const Vector3r&                shift2,
	                const bool&                    force,
	                const shared_ptr<Interaction>& c);
	virtual bool goReverse(
	        const shared_ptr<Shape>&       cm1,
	        const shared_ptr<Shape>&       cm2,
	        const State&                   state1,
	        const State&                   state2,
	        const Vector3r&                shift2,
	        const bool&                    force,
	        const shared_ptr<Interaction>& c);
	// clang-format off
	YADE_CLASS_BASE_DOC(Ig2_Facet_Sphere_ScGeom6D,Ig2_Facet_Sphere_ScGeom,"Create an interaction geometry :yref:`ScGeom6D` from :yref:`Facet` and :yref:`Sphere`, representing the Facet with a projected virtual sphere of same radius.")
	// clang-format on
	FUNCTOR2D(Facet, Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Facet, Sphere);
};
REGISTER_SERIALIZABLE(Ig2_Facet_Sphere_ScGeom6D);

class Ig2_Wall_Sphere_ScGeom : public IGeomFunctor {
public:
	virtual bool
	go(const shared_ptr<Shape>&       cm1,
	   const shared_ptr<Shape>&       cm2,
	   const State&                   state1,
	   const State&                   state2,
	   const Vector3r&                shift2,
	   const bool&                    force,
	   const shared_ptr<Interaction>& c);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Ig2_Wall_Sphere_ScGeom,IGeomFunctor,"Create/update a :yref:`ScGeom` instance representing intersection of :yref:`Wall` and :yref:`Sphere`.",
		((bool,noRatch,true,,"Avoid granular ratcheting"))
	);
	// clang-format on
	FUNCTOR2D(Wall, Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Wall, Sphere);
};
REGISTER_SERIALIZABLE(Ig2_Wall_Sphere_ScGeom);

} // namespace yade
