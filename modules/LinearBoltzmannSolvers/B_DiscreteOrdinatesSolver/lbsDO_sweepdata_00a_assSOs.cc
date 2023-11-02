#include "lbs_discrete_ordinates_solver.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_structs.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "math/Quadratures/angular_quadrature_base.h"

#include "math/Quadratures/angular_product_quadrature.h"

#define POLAR_ILLEGAL_GEOTYPE (fname + \
  ": The simulation is using polar angle aggregation for which only " \
  "certain geometry types are supported, i.e., ORTHOGONAL, DIMENSION_2 " \
  "or 3D EXTRUDED.")

#define POLAR_ONLY_PRODUCT (fname + \
  ": The simulation is using polar angle aggregation for which only " \
  "Product-type quadratures are supported.")

#define PRODUCT_QUAD_CASTING_FAILED (fname + \
  ": Casting the angular quadrature to the product quadrature base, failed.")

#define AZIMUTHAL_ILLEGAL_GEOTYPE (fname + \
  ": The simulation is using azimuthal angle aggregation for which only " \
  "ONED_SPHERICAL or TWOD_CYLINDRICAL derived geometry types are supported.")

#define AZIMUTHAL_ONLY_PRODUCT (fname + \
  ": The simulation is using azimuthal angle aggregation for which only " \
  "Product-type quadratures are supported.")

#define LogicCheck(condition, message) \
if ((condition)) \
  throw std::logic_error(fname+(message));

namespace lbs
{

//###################################################################
/**This routine groups angle-indices to groups sharing the same sweep
 * ordering. It also takes geometry into account.*/
std::pair<UniqueSOGroupings, DirIDToSOMap>
DiscreteOrdinatesSolver::
  AssociateSOsAndDirections(const chi_mesh::MeshContinuum &grid,
                            const chi_math::AngularQuadrature& quadrature,
                            const AngleAggregationType agg_type,
                            const lbs::GeometryType lbs_geo_type)
{
  const std::string fname = __FUNCTION__;

  //================================================== Checks
  LogicCheck(quadrature.omegas_.empty(),
             ": Quadrature with no omegas cannot be used.")
  LogicCheck(quadrature.weights_.empty(),
             ": Quadrature with no weights cannot be used.")

  //================================================== Build groupings
  UniqueSOGroupings unq_so_grps;
  switch (agg_type)
  {
    //=========================================== Single
    // The easiest aggregation type. Every direction
    // either has/is assumed to have a unique sweep
    // ordering. Hence there is only group holding ALL
    // the direction indices.
    case AngleAggregationType::SINGLE:
    {
      const size_t num_dirs = quadrature.omegas_.size();
      for (size_t n=0; n<num_dirs; ++n)
        unq_so_grps.push_back({n});
      break;
    }//case agg_type SINGLE

      //=========================================== Polar
      // The following conditions allow for polar
      // angle aggregation.
    case AngleAggregationType::POLAR:
    {
      //Check geometry types
      const auto grid_attribs = grid.Attributes();
      if (not(grid_attribs & chi_mesh::ORTHOGONAL or
              grid_attribs & chi_mesh::DIMENSION_2 or
              grid_attribs & chi_mesh::EXTRUDED))
        throw std::logic_error(POLAR_ILLEGAL_GEOTYPE);

      //Check quadrature type
      const auto quad_type = quadrature.type_;
      if (quad_type != chi_math::AngularQuadratureType::ProductQuadrature)
        throw std::logic_error(POLAR_ONLY_PRODUCT);

      //Process Product Quadrature
      try
      {
        typedef chi_math::ProductQuadrature ProdQuadType;
        const auto& product_quad = dynamic_cast<const ProdQuadType&>(quadrature);

        const auto num_azi = product_quad.azimu_ang_.size();
        const auto num_pol = product_quad.polar_ang_.size();

        //Make two separate list of polar angles
        //One upward-pointing and one downward
        std::vector<size_t> upward_polar_ids;
        std::vector<size_t> dnward_polar_ids;
        for (size_t p=0; p<num_pol; ++p)
          if (product_quad.polar_ang_[p] > M_PI_2)
            upward_polar_ids.push_back(p);
          else
            dnward_polar_ids.push_back(p);

        //Define lambda working for both upward and dnward polar-ids
        /**Lambda to convert indices and push it onto unq_so_grps.*/
        auto MapPolarAndAzimuthalIDs = [&product_quad, &unq_so_grps](
          const DirIDs& polar_ids,const size_t azimuthal_id)
        {
          DirIDs dir_ids;
          dir_ids.reserve(polar_ids.size());
          for (const size_t p : polar_ids)
            dir_ids.push_back(product_quad.GetAngleNum(p, azimuthal_id));
          unq_so_grps.push_back(std::move(dir_ids));
        };

        //Stack id's for all azimuthal angles
        for (size_t a=0; a<num_azi; ++a)
        {
          if (not upward_polar_ids.empty())
            MapPolarAndAzimuthalIDs(upward_polar_ids, a);
          if (not dnward_polar_ids.empty())
            MapPolarAndAzimuthalIDs(dnward_polar_ids, a);
        }//for azi-id a

      }//try product quadrature
      catch (const std::bad_cast& bc)
      {
        throw std::runtime_error(PRODUCT_QUAD_CASTING_FAILED);
      }

      break;
    }//case agg_type POLAR

      //====================================== Azimuthal
    case AngleAggregationType::AZIMUTHAL:
    {
      //Check geometry types
      if (not (lbs_geo_type == GeometryType::ONED_SPHERICAL or
               lbs_geo_type == GeometryType::TWOD_CYLINDRICAL))
        throw std::logic_error(AZIMUTHAL_ILLEGAL_GEOTYPE);

      //Check quadrature type
      const auto quad_type = quadrature.type_;
      if (quad_type != chi_math::AngularQuadratureType::ProductQuadrature)
        throw std::logic_error(AZIMUTHAL_ONLY_PRODUCT);

      //Process Product Quadrature
      try
      {
        typedef chi_math::ProductQuadrature ProdQuadType;
        const auto& product_quad = dynamic_cast<const ProdQuadType&>(quadrature);

        for (const auto& dir_set : product_quad.GetDirectionMap())
        {
          std::vector<unsigned int> group1;
          std::vector<unsigned int> group2;
          for (const auto& dir_id : dir_set.second)
            if (quadrature.abscissae_[dir_id].phi > M_PI_2)
              group1.push_back(dir_id);
            else
              group2.push_back(dir_id);

          DirIDs group1_ids(group1.begin(),group1.end());
          DirIDs group2_ids(group2.begin(),group2.end());

          unq_so_grps.push_back(std::move(group1_ids));
          unq_so_grps.push_back(std::move(group2_ids));
        }
      }//try product quadrature
      catch (const std::bad_cast& bc)
      {
        throw std::runtime_error(PRODUCT_QUAD_CASTING_FAILED);
      }

      break;
    }
    default:
      throw std::invalid_argument(fname + ": Called with UNDEFINED angle "
                                          "aggregation type.");
  }//switch angle aggregation type

  //================================================== Map directions to sweep
  //                                                   orderings
  DirIDToSOMap dir_id_to_so_map;
  {
    size_t so_grouping_id = 0;
    for (const auto& so_grouping : unq_so_grps)
    {
      for (const size_t dir_id: so_grouping)
        dir_id_to_so_map[dir_id] = so_grouping_id;

      ++so_grouping_id;
    }//for so_grouping
  }//map scope

  return {unq_so_grps, dir_id_to_so_map};
}


}//namespace lbs

