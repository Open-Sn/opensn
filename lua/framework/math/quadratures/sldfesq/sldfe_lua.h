#pragma once

/** Creates a Simplified Linear Discontinuous Finite Element (SLDFE)
 * quadrature based on Spherical Quadrilaterals (SQ). Hence SLDFE-SQ.
 * \param initial_refinement_level int Initial refinement level, \f$n\f$ to
 *        be used. The total number of angles will be \f$ 8{\times}12(n+1)^2 \f$.
 *
 * ##_
 *
 * ###Example:
 * Example with refinement level 2.
 * \code
 * pquad = CreateSLDFESQAngularQuadrature(2)
 * \endcode
 *
 * \image html "SLDFESQBasen2.png" width=500px
 * With direction points:
 * \image html "SLDFESQBasen2Oct2.png" width=500px
 *
 * \ingroup LuaSLDFESQ
 * \author Jan
 */
int CreateSLDFESQAngularQuadrature(lua_State* L);

/** Applies a local refinement of angles.
 * \param handle int. Handle to the reference quadrature.
 * \param reference_direction vec3 Reference vector. \f$ \vec{r} \f$
 * \param cone_size double Cone size in radians. \f$ \theta \f$
 * \param invert_logic bool Optional[Default:false]. If supplied, interprets
 *  SQ-splitting as when \f$|\omega \cdot \vec{r}| < \sin(\theta) \f$. Otherwise,
 *  SQs will be split if \f$ \omega \cdot \vec{r} > \cos(\theta)\f$
 *
 * ##_
 *
 * ###Example:
 * Example with refinement level 2 and a triple directional refinement:
 * \code
 * pquad = CreateSLDFESQAngularQuadrature(2)
 * chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},45.0*math.pi/180,false)
 * chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},23.0*math.pi/180,false)
 * chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},12.0*math.pi/180,false)
 * \endcode
 *
 * \image html "SLDFESQr.png" width=500px
 *
 * Example with refinement level 2 and a triple planar refinement:
 * \code
 * pquad = CreateSLDFESQAngularQuadrature(2)
 * chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},22.50*math.pi/180,true)
 * chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},11.75*math.pi/180,true)
 * chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},5.000*math.pi/180,true)
 * \endcode
 *
 * \image html "SLDFESQp.png" width=500px
 *
 * \ingroup LuaSLDFESQ
 * \author Jan
 */
int chiLocallyRefineSLDFESQAngularQuadrature(lua_State* L);

/** Outputs the quadrature information to python format.
 * \param handle int Handle to the reference quadrature.
 * \param file_name_prefix string Prefix to be used in front of file.
 *
 * ##_
 *
 * ###Example:
 * Example of printing a quadrature:
 * Example with refinement level 2 and a triple directional refinement:
 * \code
 * pquad = CreateSLDFESQAngularQuadrature(2)
 * chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},45.0*math.pi/180,false)
 * chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},23.0*math.pi/180,false)
 * chiLocallyRefineSLDFESQAngularQuadrature(pquad,{1,0,0},12.0*math.pi/180,false)
 * chiPrintToPythonSLDFESQAngularQuadrature(pquad,"YQuad_");
 * \endcode
 *
 * \image html "SLDFESQr.png" width=500px
 *
 * \ingroup LuaSLDFESQ
 * \author Jan
 */
int chiPrintToPythonSLDFESQAngularQuadrature(lua_State* L);
