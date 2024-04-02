#pragma once

#include "mpi.h"
#include "Enums.h"
#include "Error.h"

namespace mpicpp_lite {

/// Wrapper around `MPI_Group`
class Group {
public:
    enum ComparisonResult { IDENTICAL = MPI_IDENT, SIMILAR = MPI_SIMILAR, UNEQUAL = MPI_UNEQUAL };

    /// Create an empty group
    Group();

    /// Create group from an `MPI_Group`
    ///
    /// @param group `MPI_Group` used to initialize this object
    Group(const MPI_Group & group);

    /// Returns the rank of this process in the given group
    ///
    /// @return Rank of the calling process in group, or MPI_UNDEFINED if the process is not a
    ///         member
    int rank() const;

    /// Returns the size of a group
    ///
    /// @return Number of processes in the group
    int size() const;

    /// Produces a group by reordering an existing group and taking only listed members
    ///
    /// @param ranks Ranks of processes in group to appear in the new group
    /// @return New group derived from this group, in the order defined by `ranks`
    Group include(const std::vector<int> & ranks) const;

    /// Produces a group by reordering an existing group and taking only unlisted members
    ///
    /// @param ranks Array of integer ranks of processes in group not to appear in the new group
    /// @return New group derived from this group, preserving the order defined by group
    Group exclude(const std::vector<int> & ranks) const;

    /// Frees the group
    void free();

    /// Translates the rank of a process in one group to the one in another group
    ///
    /// @param in_rank Rank in this group
    /// @param out_group Destination group
    /// @return Rank in destination group, or `UNDEFINED` when no correspondence exists
    int translate_rank(int in_rank, const Group & out_group) const;

    /// Translates the ranks of processes in one group to those in another group
    ///
    /// @param in_ranks Array of valid ranks in this group
    /// @param out_group Destination group
    /// @return Array of corresponding ranks in destination group, `UNDEFINED` when no
    ///         correspondence exists
    std::vector<int> translate_ranks(const std::vector<int> & in_ranks,
                                     const Group & out_group) const;

    /// Type cast operator so we can pass this class directly into MPI calls
    operator const MPI_Group &() const { return this->group; }

public:
    /// Compares two groups
    ///
    /// @param g1 First group
    /// @param g2 Second group
    /// @return `IDENTICAL` if the order and members of the two groups are the same, `SIMILAR` if
    ///         only the members are the same, and `UNEQUAL` otherwise
    static ComparisonResult compare(const Group & g1, const Group & g2);

    /// Produces a group by combining two groups
    ///
    /// @param g1 First group
    /// @param g2 Second group
    /// @return Union group
    static Group join(const Group & g1, const Group & g2);

    /// Produces a group as the intersection of two existing groups
    ///
    /// @param g1 First group
    /// @param g2 Second group
    /// @return Intersection group
    static Group intersection(const Group & g1, const Group & g2);

    /// Makes a group from the difference of two groups
    ///
    /// @param g1 First group
    /// @param g2 Second group
    /// @return Difference group
    static Group difference(const Group & g1, const Group & g2);

private:
    MPI_Group group;
};

inline Group::Group() {}

inline Group::Group(const MPI_Group & group) : group(group) {}

inline Group
Group::include(const std::vector<int> & ranks) const
{
    MPI_Group new_group;
    MPI_CHECK(MPI_Group_incl(this->group, ranks.size(), ranks.data(), &new_group));
    return { new_group };
}

inline Group
Group::exclude(const std::vector<int> & ranks) const
{
    MPI_Group new_group;
    MPI_CHECK(MPI_Group_excl(this->group, ranks.size(), ranks.data(), &new_group));
    return { new_group };
}

inline void
Group::free()
{
    MPI_CHECK(MPI_Group_free(&this->group));
}

inline int
Group::rank() const
{
    int r;
    MPI_CHECK(MPI_Group_rank(this->group, &r));
    return r;
}

inline int
Group::size() const
{
    int sz;
    MPI_CHECK(MPI_Group_size(this->group, &sz));
    return sz;
}

inline int
Group::translate_rank(int in_rank, const Group & out_group) const
{
    int out_rank;
    MPI_CHECK(MPI_Group_translate_ranks(this->group, 1, &in_rank, out_group, &out_rank));
    return out_rank;
}

inline std::vector<int>
Group::translate_ranks(const std::vector<int> & in_ranks, const Group & out_group) const
{
    int n = in_ranks.size();
    std::vector<int> ranks2(n);
    MPI_CHECK(MPI_Group_translate_ranks(this->group, n, in_ranks.data(), out_group, ranks2.data()));
    return ranks2;
}

inline Group::ComparisonResult
Group::compare(const Group & g1, const Group & g2)
{
    int result;
    MPI_CHECK(MPI_Group_compare(g1, g2, &result));
    return static_cast<ComparisonResult>(result);
}

inline Group
Group::join(const Group & g1, const Group & g2)
{
    MPI_Group new_group;
    MPI_CHECK(MPI_Group_union(g1, g2, &new_group));
    return { new_group };
}

inline Group
Group::intersection(const Group & g1, const Group & g2)
{
    MPI_Group new_group;
    MPI_CHECK(MPI_Group_intersection(g1, g2, &new_group));
    return { new_group };
}

inline Group
Group::difference(const Group & g1, const Group & g2)
{
    MPI_Group new_group;
    MPI_CHECK(MPI_Group_difference(g1, g2, &new_group));
    return { new_group };
}

} // namespace mpicpp_lite
