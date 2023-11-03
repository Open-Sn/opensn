#pragma once

#include <vector>
#include <map>
#include <memory>
#include <array>

// ################################################################### Fwd decls
namespace chi_mesh::sweep_management
{
class SweepBoundary;
}
namespace chi_physics
{
class MultiGroupXS;
}

namespace lbs::acceleration
{

/**Boundary condition type. We essentially only support two
 * types: Dirichlet and Reflecting, the latter is covered under
 * the ROBIN-type boundary condition.*/
enum class BCType
{
  DIRICHLET = 1,
  ROBIN = 2
};

/**Simple data structure to specify boundary conditions. Its stores the
 * BC-type in `type` and an array of 3 values in `values`. For a
 * Dirichlet-BC only `values[0]` is used to specify the value of the BC.
 * For a robin boundary condition we use all 3 values in the form
\f[
a\phi + b \mathbf{n} \frac{\partial \phi}{\partial \mathbf{x}} = f
\f]
where \f$ a \f$, \f$ b \f$ and \f$ f \f$ map to `values[0]`, `values[1]` and
`values[2]`, respectively.
*/
struct BoundaryCondition
{
  BCType type = BCType::DIRICHLET;
  std::array<double, 3> values = {0, 0, 0};
};

/**For both WGDSA and TGDSA use, we define this
 * simplified data structure to hold multigroup diffusion coefficients and
 * removal cross sections only for the relevant groups. E.g., for WGDSA it
 * will only hold the cross sections for the groupset, and for TGDSA there will
 * only be one group (All groups collapsed into 1).*/
struct Multigroup_D_and_sigR
{
  std::vector<double> Dg;
  std::vector<double> sigR;
};

enum class EnergyCollapseScheme
{
  JFULL = 1,    ///< Jacobi with full conv. of within-group scattering
  JPARTIAL = 2, ///< Jacobi with partially conv. of within-group scattering
};

struct TwoGridCollapsedInfo
{
  double collapsed_D = 0.0;
  double collapsed_sig_a = 0.0;
  std::vector<double> spectrum;
};

TwoGridCollapsedInfo MakeTwoGridCollapsedInfo(const chi_physics::MultiGroupXS& xs,
                                              EnergyCollapseScheme scheme);

typedef std::shared_ptr<chi_mesh::sweep_management::SweepBoundary> SwpBndryPtr;

/**Translates sweep boundary conditions to that used in diffusion acceleration
 * methods.*/
std::map<uint64_t, BoundaryCondition>
TranslateBCs(const std::map<uint64_t, SwpBndryPtr>& sweep_boundaries,
             bool vaccum_bcs_are_dirichlet = true);

typedef std::shared_ptr<chi_physics::MultiGroupXS> MGXSPtr;

/**Makes a packaged set of XSs, suitable for diffusion, for a particular
 * set of groups.*/
std::map<int, Multigroup_D_and_sigR> PackGroupsetXS(const std::map<int, MGXSPtr>& matid_to_xs_map,
                                                    int first_grp_index,
                                                    int last_group_index);

} // namespace lbs::acceleration
