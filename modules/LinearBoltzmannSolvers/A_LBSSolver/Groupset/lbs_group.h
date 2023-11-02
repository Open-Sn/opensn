#ifndef LBS_GROUP
#define LBS_GROUP

//################################################################### Class def
namespace lbs
{
/**Object holding a grouping.*/
class LBSGroup
{
public:
  int id_;

public:
  LBSGroup() : id_(-1) {}
  explicit LBSGroup(int id) : id_(id) {}
};
} // namespace lbs

#endif