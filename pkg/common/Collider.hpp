/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<core/Bound.hpp>
#include<core/Interaction.hpp>
#include<core/GlobalEngine.hpp>

#include<pkg/common/Dispatching.hpp>


class Collider: public GlobalEngine {
	public:
		static int avoidSelfInteractionMask;
		/*! Probe the Bound on a bodies presense. Returns list of body ids with which there is potential overlap. */
		virtual  vector<Body::id_t> probeBoundingVolume(const Bound&){throw;}
		/*! Tell whether given bodies may interact, for other than spatial reasons.
		 *
		 * Concrete collider implementations should call this function if
		 * the bodies are in potential interaction geometrically. */
		static bool mayCollide(const Body*, const Body*
		#ifdef YADE_MPI
		,Body::id_t subdomain
		#endif 
		);
		/*! Invalidate all persistent data (if the collider has any), forcing reinitialization at next run.
		The default implementation does nothing, colliders should override it if it is applicable.

		Currently used from Shop::flipCell, which changes cell information for bodies.
		*/
		virtual void invalidatePersistentData(){}

		// ctor with functors for the integrated BoundDispatcher
		virtual void pyHandleCustomCtorArgs(boost::python::tuple& t, boost::python::dict& d);
		
		int get_avoidSelfInteractionMask(){return avoidSelfInteractionMask;}
		void set_avoidSelfInteractionMask(const int &v){avoidSelfInteractionMask = v;}
		
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Collider,GlobalEngine,"Abstract class for finding spatial collisions between bodies. \n\n.. admonition:: Special constructor\n\n\tDerived colliders (unless they override ``pyHandleCustomCtorArgs``) can be given list of :yref:`BoundFunctors <BoundFunctor>` which is used to initialize the internal :yref:`boundDispatcher <Collider.boundDispatcher>` instance.",
		((shared_ptr<BoundDispatcher>,boundDispatcher,new BoundDispatcher,Attr::readonly,":yref:`BoundDispatcher` object that is used for creating :yref:`bounds <Body.bound>` on collider's request as necessary.")),
		/*ctor*/,
		.add_property("avoidSelfInteractionMask",&Collider::get_avoidSelfInteractionMask,&Collider::set_avoidSelfInteractionMask , R"""(
This mask is used to avoid the interactions inside a group of particles. To do so, the particles must have the exact same mask and that mask should have one bit in common with this :yref:`avoidSelfInteractionMask<Collider.avoidSelfInteractionMask>` as for their binary representations.

In an example with 3 mask values assigned to particles: ``0b10`` = "blue", ``0b01`` = "green", ``0b11`` = "red", and a with collider.avoidSelfInteractionMask = ``0b01``:

* Green particles would not interact with blue particles because they are not mask-compatible (do not share any "1" bit in common): ``0b01 & 0b10`` is false
* Green particles would not interact with other green particles because of avoidSelfInteractionMask: ``0b01 & 0b01`` is true
* Red particles would not interact with themselves either: ``0b11 & 0b01`` is true
* For the same reason blue particles would interact with themselves because for them avoidSelfInteractionMask: ``0b10 and 0b01`` is false

For the complete list of possibilities, see :ysrc:`pkg/common/Collider.cpp#L15`.
	)""")
	);
};
REGISTER_SERIALIZABLE(Collider);
