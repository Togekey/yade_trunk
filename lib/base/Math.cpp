#include <lib/base/Math.hpp>
#include <boost/math/constants/constants.hpp>

namespace yade { // Cannot have #include directive inside.

// TODO: replace with numeric_limits.
template <> const Real Math<Real>::EPSILON        = DBL_EPSILON;
template <> const Real Math<Real>::ZERO_TOLERANCE = 1e-20;
template <> const Real Math<Real>::MAX_REAL       = DBL_MAX;

template <> const Real Math<Real>::PI             = boost::math::constants::pi<Real>();
template <> const Real Math<Real>::TWO_PI         = boost::math::constants::two_pi<Real>();
template <> const Real Math<Real>::HALF_PI        = boost::math::constants::half_pi<Real>();
template <> const Real Math<Real>::DEG_TO_RAD     = Math<Real>::PI / Real(180);
template <> const Real Math<Real>::RAD_TO_DEG     = Real(180) / Math<Real>::PI;

template<> int ZeroInitializer<int>(){ return (int)0; }
template<> Real ZeroInitializer<Real>(){ return (Real)0; }

#ifdef YADE_MASK_ARBITRARY
bool operator==(const mask_t& g, int i) { return g == mask_t(i); }
bool operator==(int i, const mask_t& g) { return g == i; }
bool operator!=(const mask_t& g, int i) { return !(g == i); }
bool operator!=(int i, const mask_t& g) { return g != i; }
mask_t operator&(const mask_t& g, int i) { return g & mask_t(i); }
mask_t operator&(int i, const mask_t& g) { return g & i; }
mask_t operator|(const mask_t& g, int i) { return g | mask_t(i); }
mask_t operator|(int i, const mask_t& g) { return g | i; }
bool operator||(const mask_t& g, bool b) { return (g != 0) || b; }
bool operator||(bool b, const mask_t& g) { return g || b; }
bool operator&&(const mask_t& g, bool b) { return (g != 0) && b; }
bool operator&&(bool b, const mask_t& g) { return g && b; }
#endif

}; // namespace yade

