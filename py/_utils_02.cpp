// This file was extracted from _utils.cpp to speed up compilation.
#include <py/_utils.hpp>
// https://codeyarns.com/2014/03/11/how-to-selectively-ignore-a-gcc-warning/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunused-function"
// Code that generates this warning, Note: we cannot do this trick in yade. If we have a warning in yade, we have to fix it! See also https://gitlab.com/yade-dev/trunk/merge_requests/73
// This method will work once g++ bug https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53431#c34 is fixed.
#include <numpy/arrayobject.h>
#pragma GCC diagnostic pop

namespace yade { // Cannot have #include directive inside.

#ifdef YADE_MPI

void initMPI()
{
	int threads;
	int rank;
	int commSize;
	MPI_Init_thread(0, 0, MPI_THREAD_SINGLE, &threads);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::cout << "myrank = " << rank << std::endl;
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	std::cout << "commSize = " << commSize << endl;
	int      color = 2; //Foam uses 1
	MPI_Comm yadeComm;  // dummy communicator;
	MPI_Comm_split(MPI_COMM_WORLD, color, rank, &yadeComm);
}

#else

void initMPI() { return; }

#endif

using math::max;
using math::min; // using inside .cpp file is ok.

py::tuple negPosExtremeIds(int axis, Real distFactor)
{
	py::tuple extrema  = Shop::aabbExtrema();
	Real      minCoord = py::extract<Real>(extrema[0][axis])(), maxCoord = py::extract<Real>(extrema[1][axis])();
	py::list  minIds, maxIds;
	for (const auto& b : *Omega::instance().getScene()->bodies) {
		shared_ptr<Sphere> s = YADE_PTR_DYN_CAST<Sphere>(b->shape);
		if (!s)
			continue;
		if (b->state->pos[axis] - s->radius * distFactor <= minCoord)
			minIds.append(b->getId());
		if (b->state->pos[axis] + s->radius * distFactor >= maxCoord)
			maxIds.append(b->getId());
	}
	return py::make_tuple(minIds, maxIds);
}

py::tuple coordsAndDisplacements(int axis, py::tuple Aabb)
{
	Vector3r bbMin(Vector3r::Zero()), bbMax(Vector3r::Zero());
	bool     useBB = py::len(Aabb) > 0;
	if (useBB) {
		bbMin = py::extract<Vector3r>(Aabb[0])();
		bbMax = py::extract<Vector3r>(Aabb[1])();
	}
	py::list retCoord, retDispl;
	FOREACH(const shared_ptr<Body>& b, *Omega::instance().getScene()->bodies)
	{
		if (useBB && !Shop::isInBB(b->state->pos, bbMin, bbMax))
			continue;
		retCoord.append(b->state->pos[axis]);
		retDispl.append(b->state->pos[axis] - b->state->refPos[axis]);
	}
	return py::make_tuple(retCoord, retDispl);
}

void setRefSe3()
{
	Scene* scene = Omega::instance().getScene().get();
	for (const auto& b : *scene->bodies) {
		b->state->refPos = b->state->pos;
		b->state->refOri = b->state->ori;
	}
	if (scene->isPeriodic) {
		scene->cell->refHSize = scene->cell->hSize;
	}
}

Real PWaveTimeStep() { return Shop::PWaveTimeStep(); };
Real RayleighWaveTimeStep() { return Shop::RayleighWaveTimeStep(); };

py::tuple interactionAnglesHistogram(int axis, int mask, size_t bins, py::tuple aabb, bool sphSph, Real minProjLen)
{
	if (axis < 0 || axis > 2)
		throw invalid_argument("Axis must be from {0,1,2}=x,y,z.");
	Vector3r bbMin(Vector3r::Zero()), bbMax(Vector3r::Zero());
	bool     useBB = py::len(aabb) > 0;
	if (useBB) {
		bbMin = py::extract<Vector3r>(aabb[0])();
		bbMax = py::extract<Vector3r>(aabb[1])();
	}
	Real              binStep = Mathr::PI / bins;
	int               axis2 = (axis + 1) % 3, axis3 = (axis + 2) % 3;
	vector<Real>      cummProj(bins, 0.);
	shared_ptr<Scene> rb = Omega::instance().getScene();
	FOREACH(const shared_ptr<Interaction>& i, *rb->interactions)
	{
		if (!i->isReal())
			continue;
		const shared_ptr<Body>&b1 = Body::byId(i->getId1(), rb), b2 = Body::byId(i->getId2(), rb);
		if (!b1->maskOk(mask) || !b2->maskOk(mask))
			continue;
		if (useBB && !Shop::isInBB(b1->state->pos, bbMin, bbMax) && !Shop::isInBB(b2->state->pos, bbMin, bbMax))
			continue;
		if (sphSph && (!dynamic_cast<Sphere*>(b1->shape.get()) || !dynamic_cast<Sphere*>(b2->shape.get())))
			continue;
		GenericSpheresContact* geom = dynamic_cast<GenericSpheresContact*>(i->geom.get());
		if (!geom)
			continue;
		Vector3r n(geom->normal);
		n[axis]   = 0.;
		Real nLen = n.norm();
		if (nLen < minProjLen)
			continue; // this interaction is (almost) exactly parallel to our axis; skip that one
		Real theta = acos(n[axis2] / nLen) * (n[axis3] > 0 ? 1 : -1);
		if (theta < 0)
			theta += Mathr::PI;
		int binNo = int(math::round(theta / binStep));
		cummProj[binNo] += nLen;
	}
	py::list val, binMid;
	for (size_t i = 0; i < (size_t)bins; i++) {
		val.append(cummProj[i]);
		binMid.append(i * binStep);
	}
	return py::make_tuple(binMid, val);
}

py::tuple bodyNumInteractionsHistogram(py::tuple aabb)
{
	Vector3r bbMin(Vector3r::Zero()), bbMax(Vector3r::Zero());
	bool     useBB = py::len(aabb) > 0;
	if (useBB) {
		bbMin = py::extract<Vector3r>(aabb[0])();
		bbMax = py::extract<Vector3r>(aabb[1])();
	}
	const shared_ptr<Scene>& rb = Omega::instance().getScene();
	vector<int>              bodyNumIntr;
	bodyNumIntr.resize(rb->bodies->size(), 0);
	int maxIntr = 0;
	FOREACH(const shared_ptr<Interaction>& i, *rb->interactions)
	{
		if (!i->isReal())
			continue;
		const Body::id_t       id1 = i->getId1(), id2 = i->getId2();
		const shared_ptr<Body>&b1 = Body::byId(id1, rb), b2 = Body::byId(id2, rb);
		if ((useBB && Shop::isInBB(b1->state->pos, bbMin, bbMax)) || !useBB) {
			if (b1->isClumpMember())
				bodyNumIntr[b1->clumpId] += 1; //count bodyNumIntr for the clump, not for the member
			else
				bodyNumIntr[id1] += 1;
		}
		if ((useBB && Shop::isInBB(b2->state->pos, bbMin, bbMax)) || !useBB) {
			if (b2->isClumpMember())
				bodyNumIntr[b2->clumpId] += 1; //count bodyNumIntr for the clump, not for the member
			else
				bodyNumIntr[id2] += 1;
		}
		maxIntr = max(max(maxIntr, bodyNumIntr[b1->getId()]), bodyNumIntr[b2->getId()]);
		if (b1->isClumpMember())
			maxIntr = max(maxIntr, bodyNumIntr[b1->clumpId]);
		if (b2->isClumpMember())
			maxIntr = max(maxIntr, bodyNumIntr[b2->clumpId]);
	}
	vector<int> bins;
	bins.resize(maxIntr + 1, 0);
	for (size_t id = 0; id < bodyNumIntr.size(); id++) {
		const shared_ptr<Body>& b = Body::byId(id, rb);
		if (b) {
			if (bodyNumIntr[id] > 0)
				bins[bodyNumIntr[id]] += 1;
			// 0 is handled specially: add body to the 0 bin only if it is inside the bb requested (if applicable)
			// otherwise don't do anything, since it is outside the volume of interest
			else if (((useBB && Shop::isInBB(b->state->pos, bbMin, bbMax)) || !useBB) && !(b->isClumpMember()))
				bins[0] += 1;
		}
	}
	py::list count, num;
	for (size_t n = 0; n < bins.size(); n++) {
		if (bins[n] == 0)
			continue;
		num.append(n);
		count.append(bins[n]);
	}
	return py::make_tuple(num, count);
}

Vector3r inscribedCircleCenter(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2) { return Shop::inscribedCircleCenter(v0, v1, v2); }
py::dict getViscoelasticFromSpheresInteraction(Real tc, Real en, Real es)
{
	shared_ptr<ViscElMat> b = shared_ptr<ViscElMat>(new ViscElMat());
	Shop::getViscoelasticFromSpheresInteraction(tc, en, es, b);
	py::dict d;
	d["kn"] = b->kn;
	d["cn"] = b->cn;
	d["ks"] = b->ks;
	d["cs"] = b->cs;
	return d;
}
/* reset highlight of all bodies */
void highlightNone()
{
	for (const auto& b : *Omega::instance().getScene()->bodies) {
		if (!b->shape)
			continue;
		b->shape->highlight = false;
	}
}

/*!Sum moments acting on given bodies
 *
 * @param ids is the calculated bodies ids
 * @param axis is the direction of axis with respect to which the moment is calculated.
 * @param axisPt is a point on the axis.
 *
 * The computation is trivial: moment from force is is by definition r×F, where r
 * is position relative to axisPt; moment from moment is m; such moment per body is
 * projected onto axis.
 */
Real sumTorques(py::list ids, const Vector3r& axis, const Vector3r& axisPt)
{
	shared_ptr<Scene> rb = Omega::instance().getScene();
	rb->forces.sync();
	Real   ret = 0;
	size_t len = py::len(ids);
	for (size_t i = 0; i < len; i++) {
		const Body*     b = (*rb->bodies)[py::extract<int>(ids[i])].get();
		const Vector3r& m = rb->forces.getTorque(b->getId());
		const Vector3r& f = rb->forces.getForce(b->getId());
		Vector3r        r = b->state->pos - axisPt;
		ret += axis.dot(m + r.cross(f));
	}
	return ret;
}
/* Sum forces acting on bodies within mask.
 *
 * @param ids list of ids
 * @param direction direction in which forces are summed
 *
 */
Real sumForces(py::list ids, const Vector3r& direction)
{
	shared_ptr<Scene> rb = Omega::instance().getScene();
	rb->forces.sync();
	Real   ret = 0;
	size_t len = py::len(ids);
	for (size_t i = 0; i < len; i++) {
		Body::id_t      id = py::extract<int>(ids[i]);
		const Vector3r& f  = rb->forces.getForce(id);
		ret += direction.dot(f);
	}
	return ret;
}

/* Sum force acting on facets given by their ids in the sense of their respective normals.
   If axis is given, it will sum forces perpendicular to given axis only (not the the facet normals).
*/
Real sumFacetNormalForces(vector<Body::id_t> ids, int axis)
{
	shared_ptr<Scene> rb = Omega::instance().getScene();
	rb->forces.sync();
	Real ret = 0;
	FOREACH(const Body::id_t id, ids)
	{
		Facet* f = YADE_CAST<Facet*>(Body::byId(id, rb)->shape.get());
		if (axis < 0)
			ret += rb->forces.getForce(id).dot(f->normal);
		else {
			Vector3r ff = rb->forces.getForce(id);
			ff[axis]    = 0;
			ret += ff.dot(f->normal);
		}
	}
	return ret;
}


/* Tell us whether a point lies in polygon given by array of points.
 *  @param xy is the point that is being tested
 *  @param vertices is Numeric.array (or list or tuple) of vertices of the polygon.
 *         Every row of the array is x and y coordinate, numer of rows is >= 3 (triangle).
 *
 * Copying the algorithm from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
 * is gratefully acknowledged:
 *
 * License to Use:
 * Copyright (c) 1970-2003, Wm. Randolph Franklin
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
 *   2. Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
 *   3. The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission. 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * http://numpy.scipy.org/numpydoc/numpy-13.html told me how to use Numeric.array from c
 */
bool pointInsidePolygon(py::tuple xy, py::object vertices)
{
	Real           testx = py::extract<Real>(xy[0])(), testy = py::extract<Real>(xy[1])();
	char**         vertData;
	int            rows, cols;
	PyArrayObject* vert   = (PyArrayObject*)vertices.ptr();
	int            result = PyArray_As2D((PyObject**)&vert /* is replaced */, &vertData, &rows, &cols, PyArray_DOUBLE);
	if (result != 0)
		throw invalid_argument("Unable to cast vertices to 2d array");
	if (cols != 2 || rows < 3)
		throw invalid_argument("Vertices must have 2 columns (x and y) and at least 3 rows.");
	int  i /*current node*/, j /*previous node*/;
	bool inside = false;
	for (i = 0, j = rows - 1; i < rows; j = i++) {
		double vx_i = *(double*)(vert->data + i * vert->strides[0]), vy_i = *(double*)(vert->data + i * vert->strides[0] + vert->strides[1]),
		       vx_j = *(double*)(vert->data + j * vert->strides[0]), vy_j = *(double*)(vert->data + j * vert->strides[0] + vert->strides[1]);
		if (((vy_i > testy) != (vy_j > testy)) && (testx < (vx_j - vx_i) * (testy - vy_i) / (vy_j - vy_i) + vx_i))
			inside = !inside;
	}
	Py_DECREF(vert);
	return inside;
}

/*! Compute area of a simple 2d polygon, using the Surveyor's formula.
	http://en.wikipedia.org/wiki/Polygon

	The first and last points shouldn't be given twice.
*/
Real simplePolygonArea2d(vector<Vector2r> P)
{
	Real   ret = 0.;
	size_t n   = P.size();
	for (size_t i = 0; i < n - 1; i++) {
		ret += P[i][0] * P[i + 1][1] - P[i + 1][0] * P[i][1];
	}
	ret += P[n - 1][0] * P[0][1] - P[0][0] * P[n - 1][1];
	return math::abs(ret / 2.);
}

/* Compute area of convex hull when when taking (swept) spheres crossing the plane at coord, perpendicular to axis.

	All spheres that touch the plane are projected as hexagons on their circumference to the plane.
	Convex hull from this cloud is computed.
	The area of the hull is returned.

*/
Real approxSectionArea(Real coord, int axis)
{
	std::list<Vector2r> cloud;
	if (axis < 0 || axis > 2)
		throw invalid_argument("Axis must be ∈ {0,1,2}");
	const int  ax1 = (axis + 1) % 3, ax2 = (axis + 2) % 3;
	const Real sqrt3 = sqrt(3);
	Vector2r   mm, mx;
	int        i = 0;
	for (const auto& b : *Omega::instance().getScene()->bodies) {
		Sphere* s = dynamic_cast<Sphere*>(b->shape.get());
		if (!s)
			continue;
		const Vector3r& pos(b->state->pos);
		const Real      r(s->radius);
		if ((pos[axis] > coord && (pos[axis] - r) > coord) || (pos[axis] < coord && (pos[axis] + r) < coord))
			continue;
		Vector2r c(pos[ax1], pos[ax2]);
		cloud.push_back(c + Vector2r(r, 0.));
		cloud.push_back(c + Vector2r(-r, 0.));
		cloud.push_back(c + Vector2r(r / 2., sqrt3 * r));
		cloud.push_back(c + Vector2r(r / 2., -sqrt3 * r));
		cloud.push_back(c + Vector2r(-r / 2., sqrt3 * r));
		cloud.push_back(c + Vector2r(-r / 2., -sqrt3 * r));
		if (i++ == 0) {
			mm = c, mx = c;
		}
		mm = Vector2r(min(c[0] - r, mm[0]), min(c[1] - r, mm[1]));
		mx = Vector2r(max(c[0] + r, mx[0]), max(c[1] + r, mx[1]));
	}
	if (cloud.size() == 0)
		return 0;
	ConvexHull2d     ch2d(cloud);
	vector<Vector2r> hull = ch2d();
	return simplePolygonArea2d(hull);
}
/* Find all interactions deriving from NormShearPhys that cross plane given by a point and normal
	(the normal may not be normalized in this case, though) and sum forces (both normal and shear) on them.

	Returns a 3-tuple with the components along global x,y,z axes, which can be viewed as "action from lower part, towards
	upper part" (lower and upper parts with respect to the plane's normal).

	(This could be easily extended to return sum of only normal forces or only of shear forces.)
*/
Vector3r forcesOnPlane(const Vector3r& planePt, const Vector3r& normal)
{
	Vector3r ret(Vector3r::Zero());
	Scene*   scene = Omega::instance().getScene().get();
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal())
			continue;
		NormShearPhys* nsi = dynamic_cast<NormShearPhys*>(I->phys.get());
		if (!nsi)
			continue;
		Vector3r pos1, pos2;
		pos1      = Body::byId(I->getId1(), scene)->state->pos;
		pos2      = Body::byId(I->getId2(), scene)->state->pos;
		Real dot1 = (pos1 - planePt).dot(normal), dot2 = (pos2 - planePt).dot(normal);
		if (dot1 * dot2 > 0)
			continue; // both (centers of) bodies are on the same side of the plane=> this interaction has to be disregarded
		// if pt1 is on the negative plane side, d3dg->normal.Dot(normal)>0, the force is well oriented;
		// otherwise, reverse its contribution. So that we return finally
		// Sum [ ( normal(plane) dot normal(interaction= from 1 to 2) ) "nsi->force" ]
		ret += (dot1 < 0. ? 1. : -1.) * (nsi->normalForce + nsi->shearForce);
	}
	return ret;
}

} // namespace yade

