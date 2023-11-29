#pragma once

/** Creates an angular quadrature.
 *
 * \param azimuthal_angles array A lua table with N entries each being an azimuthal
 *                               angle.
 * \param polar_angles array A lua table with N entries each being a polar angle
 *                               angle.
 * \param weight array A lua table with N entries each being a quadrature weight.
 *
 *
 * \return Returns a unique handle to the created angular quadrature
 *
 * \ingroup LuaQuadrature
 * \author Jan
 */
int chiCreateCustomAngularQuadrature(lua_State* L);

/** Creates a curvilinear product quadrature suitable for cylindrical
 * geometries.
 *
 * \param QuadratureType int Quadrature identifier.
 * \param values varying Varying options based on the quadrature type.
 *
 * ##_
 *
 * ###QuadratureType:\n
 * GAUSS_LEGENDRE_CHEBYSHEV\n
 *   Gauss-Legendre quadrature for the polar angle and Gauss-Chebyshev quadrature
 *   for the azimuthal angle.
 *   Arguments for this quadrature type, in order:
 *   - Np : (int) number of polar angles
 *   - Na : (int) number of azimuthal angles (unique number at each polar level),
 * or (table<int>) number of azimuthal angles (diverse number at each polar level)
 *   - verbose : (bool) verbosity flag (optional).
 *
 * ###QuadratureType:\n
 * GAUSS_LEGENDRE_LEGENDRE\n
 *   Gauss-Legendre quadrature for the polar angle and Gauss-Legendre quadrature
 *   for the azimuthal angle.
 *   Arguments for this quadrature type, in order:
 *   - Np : (int) number of polar angles
 *   - Na : (int) number of azimuthal angles (unique number at each polar level),
 * or (table<int>) number of azimuthal angles (diverse number at each polar level)
 *   - verbose : (bool) verbosity flag (optional).
 *
 *
 * \return Returns a unique handle to the created product quadrature rule
 *
 * \ingroup LuaQuadrature
 */
int chiCreateCylindricalProductQuadrature(lua_State* L);
int chiCreateSphericalProductQuadrature(lua_State* L);

/** Creates a Product-quadrature.
 *
 * \param QuadratureType int Quadrature identifier.
 * \param values varying Varying options based on the quadrature type.
 *
 * ##_
 *
 * ###QuadratureType:
 * GAUSS_LEGENDRE\n
 *  Gauss-Legendre quadrature for the polar angles and no quadrature rule
 *  for the azimuthal angle. Suitable only for 1D simulations. Expects
 *  to be followed by the number of angles Np. Optionally a verbosity flag
 *  can be added.\n\n
 *
 * GAUSS_LEGENDRE_LEGENDRE\n
 *  Gauss-Legendre quadrature for both the polar and azimuthal dimension.
 *  Expects to be followed by number of Azimuthal and Polar angles.
 *  Optionally a verbosity flag can be added.\n\n
 *
 * GAUSS_LEGENDRE_CHEBYSHEV\n
 *  Gauss-Legendre quadrature for the polar angle but Gauss-Chebyshev
 *  for the azimuthal angle.
 *  Expects to be followed by number of Azimuthal and Polar angles.
 *  Optionally a verbosity flag can be added.\n\n
 *
 * CUSTOM_QUADRATURE\n
 *  Expects to be followed by three lua tables. The first table is an array,
 *  of length Na, of the azimuthal angles (radians). The second table is an array,
 *  of length Np, of the polar angles (radians). The third table is an array, of
 *  length Na*Np, and contains the weight associated with each angle pair.
 *  Optionally a verbosity flag can be added.\n\n
 *
 *
 * \return Returns a unique handle to the created product quadrature rule
 *
 * \ingroup LuaQuadrature
 * \author Jan
 */
int CreateProductQuadrature(lua_State* L);

/** Creates a quadrature.
 *
 * \param QuadratureType int Quadrature identifier.
 * \param NumberOfPoints int Number of quadrature points.
 * \param VerboseFlag bool As the name implies. Default: false.
 *
 * ##_
 *
 * ###QuadratureType:\n
 *  GAUSS_LEGENDRE = Gauss-Legendre quadrature.
 *  GAUSS_CHEBYSHEV = Gauss-Chebyshev quadrature.
 *
 * \return Returns a unique handle to the created quadrature rule
 *
 * \ingroup LuaQuadrature
 * \author Jan
 */
int chiCreateLineQuadrature(lua_State* L);

/** Get the values of a product quadrature
 *
 * \param QuadHandle int Handle to an existing product quadrature.
 *
 * \return Table A lua table with each entry being another table with entries
 *        .weight .polar .azimuthal.
 *
 * \ingroup LuaQuadrature
 * \author Jan
 */
int chiGetProductQuadrature(lua_State* L);

/**Optimizes the indicated angular quadrature for polar symmetry.
 *
 * \param handle        int. Handle to the quadrature to be optimized.
 * \param normalization double. (Optional) The normalization to be applied to the
 *                      modified quadrature. Any negative number will inhibit
 *                      renormalization. [Default=-1.0]
 *
 *  ## _
 *
 *  ###Example:
 *  Example:
 * \code
 * pquad = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 1)
 * chiOptimizeAngularQuadratureForPolarSymmetry(pqaud, 4.0*math.pi)
 * \endcode
 *
 * \ingroup LuaQuadrature
 * \author Jan
 */
int chiOptimizeAngularQuadratureForPolarSymmetry(lua_State* L);
