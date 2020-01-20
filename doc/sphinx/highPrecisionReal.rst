.. _highPrecisionReal:

***************************
High precision calculations
***************************

Yade supports high precision and arbitrary precision ``Real`` type for performing calculations. All tests and checks pass but still the current support is considered experimental.
The backend library is `boost <https://github.com/boostorg/multiprecision>`__ `multiprecision <https://www.boost.org/doc/libs/1_72_0/libs/multiprecision/doc/html/index.html>`__
along with corresponding `boost <https://github.com/boostorg/math>`__ `math toolkit <https://www.boost.org/doc/libs/1_72_0/libs/math/doc/html/index.html>`__.

The supported types are following:

=============================================== =============== =============================== ==================================================================
type						bits		decimal places			notes
=============================================== =============== =============================== ==================================================================
``float``					``32``		``6``				hardware accelerated (not useful, it is only for testing purposes)
``double``					``64``		``15``				hardware accelerated
``long double``					``80``		``18``				hardware accelerated
``boost::multiprecision::float128``		``128``		``33``				depending on processor type it may be hardware accelerated
``boost::multiprecision::mpfr``			``Nbit``	``Nbit/(log(2)/log(10))``	uses external mpfr library
``boost::multiprecision::cpp_bin_float``	``Nbit``	``Nbit/(log(2)/log(10))``	uses boost only, but is slower
=============================================== =============== =============================== ==================================================================

The last two types are arbitrary precision, and their number of bits ``Nbit`` or decimal places is specified as argument during compilation.

Installation
===========================================

The precompiled yade package uses ``double`` type by default. In order to use high precision type yade has to be compiled and installed from source code by following the
regular :ref:`installation instructions <install-from-source-code>`. With extra following caveats:

1. Additional packages are required to be installed: ``python-mpmath`` ``libmpfr-dev`` ``libmpfrc++-dev`` (the ``mpfr`` related
   packages are necessary only to use ``boost::multiprecision::mpfr`` type).

2. The debian/ubuntu package ``python-minieigen`` must be of version 0.50.4 or higher. This package hasn't been prepared yet, hence currently following solutions are possible:

	1. Patch the older ``python-minieigen`` debian package with `these patches <https://gitlab.com/cosurgi/minieigen-real/tree/master/patches>`__, the patches already included in the debian package are also required.
	2. Or copy `the patched minieigen files <https://gitlab.com/cosurgi/minieigen-real/tree/master/minieigen-local>`__ into local directory and pass the extra cmake argument ``MINIEIGEN_INCLUDE_PATH=…``, relative to this path the sources should reside inside ``minieigen`` subdirectory. The path ``/usr/include/minieigen`` is used by default to find the minieigen sources.

3. A g++ compiler version 9.2.1 or higher is required. It shall be noted that upgrading compiler only on an existing linux installation is difficult and it is not recommended. A simpler solution is to upgrade current linux installation.

4. During cmake invocation specify:

	1. either number of bits as ``REAL_PRECISION_BITS=……``
	2. or number of requested decimal places as ``REAL_DECIMAL_PLACES=…``, but not both
	3. to use MPFR specify ``ENABLE_MPFR=ON``.

   The arbitrary precision (``mpfr`` or ``cpp_bin_float``) types are used only when more than 128 bits or more than 39 decimal places are requested. In such case if ``ENABLE_MPFR=OFF`` then
   the slower ``cpp_bin_float`` type is used. The difference in decimal places between 39 and 33 stems from the fact that `15 bits are used for exponent <https://en.wikipedia.org/wiki/Quadruple-precision_floating-point_format>`__.

.. FIXME : MINIEIGEN_PATH is not used by cmake currently


Supported modules
===========================================

During :ref:`compilation <yadeCompilation>` several yade modules can be enabled or disabled by passing an ``ENABLE_*`` command line argument to cmake.
The following table lists which modules are currently supported using high precision:

=========================================== ============================= ========================
``ENABLE_*`` module name                    default cmake setting         module support status
=========================================== ============================= ========================
``ENABLE_GUI``                              ``ON``                        native support [#supp1]_
``ENABLE_CGAL``                             ``ON``                        native support [#supp1]_
``ENABLE_VTK``                              ``ON``                        supported [#supp3]_
``ENABLE_OPENMP``                           ``ON if Nbits<=128``          conditionally supported [#supp5]_
``ENABLE_MPI``                              ``OFF``                       not tested [#supp6]_
``ENABLE_GTS``                              ``ON``                        supported [#supp2]_
``ENABLE_GL2PS``                            ``ON``                        supported [#supp2]_
``ENABLE_LINSOLV``                          ``OFF``                       not supported [#supp7]_
``ENABLE_PFVFLOW``                          ``OFF``                       not supported [#supp7]_
``ENABLE_TWOPHASEFLOW``                     ``OFF``                       not supported [#supp7]_
``ENABLE_THERMAL``                          ``OFF``                       not supported [#supp7]_
``ENABLE_LBMFLOW``                          ``ON``                        supported [#supp2]_
``ENABLE_SPH``                              ``OFF``                       not tested [#supp9]_
``ENABLE_LIQMIGRATION``                     ``OFF``                       not tested [#supp9]_
``ENABLE_MASK_ARBITRARY``                   ``OFF``                       not tested [#supp9]_
``ENABLE_PROFILING``                        ``OFF``                       not tested [#supp9]_
``ENABLE_POTENTIAL_BLOCKS``                 ``OFF``                       not supported [#supp8]_
``ENABLE_POTENTIAL_PARTICLES``              ``ON``                        supported [#supp4]_
``ENABLE_DEFORM``                           ``OFF``                       not tested [#supp9]_
``ENABLE_FEMLIKE``                          ``ON``                        supported [#supp2]_
``ENABLE_OAR``                              ``OFF``                       not tested [#supp9]_
``ENABLE_ASAN``                             ``OFF``                       supported [#supp2]_
``ENABLE_MPFR``                             ``OFF``                       native support [#supp1]_
=========================================== ============================= ========================

.. rubric:: Footnotes

.. [#supp1] this feature is supported natively, which means that specific numerical traits were written :ysrc:`for Eigen <lib/high-precision/EigenNumTraits.hpp>` and :ysrc:`for CGAL <lib/high-precision/CgalNumTraits.hpp>`, as well as :ysrc:`specific GUI <gui/qt5/SerializableEditor.py>` and :ysrc:`python support <lib/high-precision/ToFromPythonConverter.hpp>` were added.

.. [#supp2] the module was tested, the ``yade --test`` and ``yade --check`` pass, as well as most of examples are working. But it hasn't been tested extensively for all possible use cases.

.. [#supp3] supported via :ysrc:`compatibility layer <lib/compatibility/VTKCompatibility.hpp>` which converts all numbers down to ``double`` type. See :ref:`below <vtkRealCompatibility>`.

.. [#supp4] the module is enabled by default, the ``yade --test`` and ``yade --check`` pass, as well as most of examples are working. However the calculations are performed at lower ``double`` precision. A wrapper/converter layer for :ysrc:`LAPACK library <lib/compatibility/LapackCompatibility.hpp>` has been implemented. To make it work with full precision these routines should be reimplemented using Eigen.

.. [#supp5] experimentally supported for ``long double`` and ``float128``. Not supported for ``mpfr`` and ``cpp_bin_float``. File :ysrc:`lib/base/openmp-accu.hpp` needs further testing. If in doubt, compile yade with ``ENABLE_OPENMP=OFF``.

.. [#supp6] MPI support has not been tested and sending data over network hasn't been tested yet.

.. [#supp7] not supported, the code uses external cholmod library which supports only ``double`` type. To make it work a native Eigen solver for linear equations should be used.

.. [#supp8] potential blocks use external library coinor for linear programming, this library uses ``double`` type only. To make it work a linear programming routine has to be implemented using Eigen, or a similar converter/wrapper can be used as for :ysrc:`LAPACK library <lib/compatibility/LapackCompatibility.hpp>`.

.. [#supp9] this feature is ``OFF`` by default, the support of this feature has not been tested.

The unsupported modules are automatically disabled during the cmake stage.

Compatibility with other libraries
===========================================

Python compatibility with high precision types
----------------------------------------------

Eigen compatibility with high precision types
----------------------------------------------


CGAL compatibility with high precision types
----------------------------------------------


.. _vtkRealCompatibility:

VTK compatibility with high precision types
-------------------------------------------

Since VTK is only used to record results for later viewing inn other software, such as `paraview <https://www.paraview.org/>`__, the recording all decimal places does not seem to be necessary.
Hence all recording commands in ``C++`` convert ``Real`` type down to ``double`` using ``static_cast<double>`` command. This has been implemented via a set of macros in :ysrc:`VTKCompatibility.hpp <lib/compatibility/VTKCompatibility.hpp>` and :ysrc:`DoubleCompatibility.hpp <lib/compatibility/DoubleCompatibility.hpp>` files.

.. comment .......................................


LAPACK compatibility with high precision types
----------------------------------------------

