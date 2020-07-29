// This file was extracted from _utils.cpp to speed up compilation.

namespace yade { // Cannot have #include directive inside.

py::tuple negPosExtremeIds(int axis, Real distFactor);
py::tuple coordsAndDisplacements(int axis, py::tuple Aabb);
void      setRefSe3();
Real      PWaveTimeStep();
Real      RayleighWaveTimeStep();
py::tuple interactionAnglesHistogram(int axis, int mask, size_t bins, py::tuple aabb, bool sphSph, Real minProjLen);
py::tuple bodyNumInteractionsHistogram(py::tuple aabb);
Vector3r  inscribedCircleCenter(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2);
py::dict  getViscoelasticFromSpheresInteraction(Real tc, Real en, Real es);
void      highlightNone();
Real      sumTorques(py::list ids, const Vector3r& axis, const Vector3r& axisPt);
Real      sumForces(py::list ids, const Vector3r& direction);
Real      sumFacetNormalForces(vector<Body::id_t> ids, int axis);
bool      pointInsidePolygon(py::tuple xy, py::object vertices);
Real      simplePolygonArea2d(vector<Vector2r> P);
Real      approxSectionArea(Real coord, int axis);
Vector3r  forcesOnPlane(const Vector3r& planePt, const Vector3r& normal);

} // namespace yade

