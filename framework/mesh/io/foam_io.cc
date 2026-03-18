// SPDX-FileCopyrightText: 2025 The OpenSn Authors
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/utils.h"
#include <algorithm>
#include <fstream>

namespace opensn
{
namespace
{
// Utility functions

// Skip whitespace and both // and /* */ comments
inline void
SkipWSAndComments(std::istream& s)
{
  while (true)
  {
    s >> std::ws;
    int c = s.peek();
    if (c != '/')
      return;
    s.get(); // '/'
    int c2 = s.peek();
    if (c2 == '/')
    { // // ...
      s.get();
      s.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }
    if (c2 == '*')
    { // /* ... */
      s.get();
      int prev = 0, cur = 0;
      while ((cur = s.get()) != std::char_traits<char>::eof())
      {
        if (prev == '*' and cur == '/')
          break;
        prev = cur;
      }
      continue;
    }
    s.unget(); // not a comment; put back '/'
    return;
  }
}

// check match for expected character
inline bool
ExpectChar(std::istream& s, char expect_c)
{
  SkipWSAndComments(s);
  int c = s.get();
  return c == expect_c;
}

// points: formatted like "(x y z)" or "x y z" (works inline as well)
inline bool
ReadPoint(std::istream& in, double& x, double& y, double& z)
{
  in >> std::ws;
  int c = in.peek();
  bool had_paren = false;
  if (c == '(')
  {
    in.get();
    had_paren = true;
  }

  if (not(in >> x >> y >> z))
    return false;

  if (had_paren)
  {
    in >> std::ws;
    int d = in.get();
    if (d != ')')
      return false; // malformed "(x y z"
  }
  return true;
}

// faces : "n(v0 v1 ...)" possibly compact or multiline (comments allowed)
inline bool
ReadFace(std::istream& in, std::vector<int>& out)
{
  SkipWSAndComments(in);

  // read vertex count n
  int n = 0;
  if (not(in >> n) or n < 0)
    return false;

  // find opening '('
  if (not ExpectChar(in, '('))
    return false;

  // read n vertex labels
  out.clear();
  out.reserve(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
  {
    SkipWSAndComments(in);
    int v = 0;
    if (not(in >> v))
      return false;
    out.push_back(v);
  }

  // find closing ')'
  return ExpectChar(in, ')');
}

// used to read the cell idx in the owner and neighbor files
// current idx in list corresponds to face and the idx parsed
// is this face's cell idx
inline bool
ReadFace2CellID(std::istream& in, int& v)
{
  SkipWSAndComments(in);
  return static_cast<bool>(in >> v);
}

// Skip OpenFOAM header block until data is parsed
void
SkipFoamHeader(std::istream& in)
{
  std::string line;
  while (std::getline(in, line))
  {
    auto tmp = StringTrim(line);
    if (tmp.empty() or StartsWith(tmp, "//"))
      continue;
    if (not tmp.empty() and (std::isdigit(static_cast<unsigned char>(tmp[0])) or tmp[0] == '-'))
    {
      // rewind to start of this line
      in.seekg(-static_cast<std::streamoff>(line.size()) - 1, std::ios_base::cur);
      return;
    }
  }
  in.clear();
}

// skip chars until target char is encountered
inline bool
SkipUntil(std::istream& in, char target)
{
  int ch = 0;
  while ((ch = in.get()) != EOF)
    if (static_cast<char>(ch) == target)
      return true;
  return false;
}

// read exactly 'n' integer entries in list ( for cellZones )
inline bool
ReadNInts(std::istream& in, int n, std::vector<int>& out)
{
  out.reserve(out.size() + static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
  {
    int v = 0;
    if (not(in >> v))
      return false;
    out.emplace_back(v);
  }
  return true;
}

// containers
struct CellZoneEntry
{
  std::string name;
  std::vector<int> cells;
};

struct BoundaryPatchEntry
{
  std::string name;
  std::string type;
  int n_faces = 0;
  int start_face = 0;
};

struct FaceLocation
{
  size_t cell_id = 0;
  size_t local_face_id = 0;
  bool valid = false;
};

// generalized function to parse/select OpenFOAM Lists (data structure in file)
// Supports both paren lists:  N ( item0 item1 ... )
// and uniform brace lists:    N { value }   ==> repeats 'value' N times.
template <class ItemReader>
void
ReadFoamList(std::istream& in,
             size_t& count_out,
             ItemReader item_fn,
             const std::string& fname,
             const char* name)
{

  // read count
  SkipWSAndComments(in);
  if (not(in >> count_out))
    throw std::logic_error(fname + ": Bad " + std::string(name) + " count.");

  // detect opener: '(' ( list ) or '{' ( uniform repeat )
  SkipWSAndComments(in);
  int opener = in.peek();
  if (opener != '(' and opener != '{')
    throw std::logic_error(fname + ": Expected '(' or '{' after " + std::string(name) + " count.");
  in.get(); // consume opener

  const bool brace_uniform = (opener == '{');

  if (not brace_uniform)
  {
    // Parenthesis list: read count_out items normally
    for (size_t i = 0; i < count_out; ++i)
      if (not item_fn(in, i))
        throw std::logic_error(fname + ": Failed reading " + std::string(name) + " item " +
                               std::to_string(i));

    // find closing ')'
    SkipWSAndComments(in);
    if (in.peek() == ')')
      in.get();
    else
    {
      // catch to parse for ')' if earlier statement failed
      // ( i.e. due to comments )
      bool closed = false;
      while (in)
      {
        int ch = in.get();
        if (not in)
          break;
        if (ch == ')')
        {
          closed = true;
          break;
        }
      }
      if (not closed)
        throw std::logic_error(fname + ": Missing ')' for " + std::string(name));
    }
  }
  else
  {
    // Brace uniform list: read the sole entry, then repeat it count_out times
    // Remember where the value starts to re-read it each iteration.
    SkipWSAndComments(in);
    std::streampos value_pos = in.tellg();

    // We must call item_fn exactly count_out times. To do that, we rewind
    // to the start of the value before each call so it reads the same token again.
    for (size_t i = 0; i < count_out; ++i)
    {
      in.clear();
      in.seekg(value_pos);
      if (not item_fn(in, i))
        throw std::logic_error(fname + ": Failed reading repeated " + std::string(name) +
                               " value (brace syntax) at item " + std::to_string(i));
    }

    // catch closure char ' } '
    SkipWSAndComments(in);
    if (in.peek() == '}')
      in.get();
    else
    {
      // catch for ' } ' in event that comments or extra whitespace exists
      bool closed = false;
      while (in)
      {
        int ch = in.get();
        if (not in)
          break;
        if (ch == '}')
        {
          closed = true;
          break;
        }
      }
      if (not closed)
        throw std::logic_error(fname + ": Missing '}' for " + std::string(name));
    }
  }
}

// Read next non-empty, non-comment line
inline bool
ReadNextContentLine(std::istream& s, std::string& out)
{
  std::string l;
  while (std::getline(s, l))
  {
    l = StringTrim(l);
    if (l.empty() or StartsWith(l, "//"))
      continue;
    out = std::move(l);
    return true;
  }
  return false;
}

// Extract first integer token from a line
inline bool
ExtractFirstIntToken(const std::string& line, int& value)
{
  std::istringstream iss(line);
  std::string tok;
  while (iss >> tok)
  {
    while (not tok.empty() and (tok.back() == '(' or tok.back() == ')' or tok.back() == ';'))
      tok.pop_back();

    if (tok == "List<label>" or tok == "labelList" or tok == "label")
      continue;

    char* endp = nullptr;
    long val = std::strtol(tok.c_str(), &endp, 10);
    if (endp and *endp == '\0')
    {
      value = static_cast<int>(val);
      return true;
    }
  }
  return false;
}

// Parse "keyword value" for patch type in OpenFOAM's boundary file
inline bool
ReadBoundaryKeywordValue(const std::string& line, std::string& key, std::string& value)
{
  std::string tmp = StringTrim(line);
  if (tmp.empty() or StartsWith(tmp, "//"))
    return false;

  if (not tmp.empty() and tmp.back() == ';')
    tmp.pop_back();

  std::istringstream iss(tmp);
  if (not(iss >> key))
    return false;

  std::getline(iss, value);
  value = StringTrim(value);

  return not key.empty();
}

std::vector<BoundaryPatchEntry>
ReadBoundaryFile(const std::filesystem::path& path, const std::string& fname)
{
  std::ifstream in(path);
  if (not in.is_open())
    return {}; // no boundary file is allowed

  SkipFoamHeader(in);

  int npatches_i = 0;
  {
    std::string line;
    if (not ReadNextContentLine(in, line) or not ExtractFirstIntToken(line, npatches_i))
      throw std::logic_error(fname + ": Bad or missing boundary patch count.");

    if (not SkipUntil(in, '('))
      throw std::logic_error(fname + ": Missing '(' starting boundary patch list.");
  }

  const auto npatches = static_cast<size_t>(npatches_i);

  std::vector<BoundaryPatchEntry> patches;
  patches.reserve(npatches);

  for (size_t p = 0; p < npatches; ++p)
  {
    std::string line;
    if (not ReadNextContentLine(in, line))
      throw std::logic_error(fname + ": Unexpected EOF while reading boundary patch #" +
                             std::to_string(p));

    if (line == ")")
      throw std::logic_error(fname + ": Unexpected ')' while reading boundary patch #" +
                             std::to_string(p));

    BoundaryPatchEntry patch;
    patch.name = line;

    if (not ReadNextContentLine(in, line) or line != "{")
      throw std::logic_error(fname + ": Expected '{' after boundary patch name '" + patch.name +
                             "'.");

    bool found_n_faces = false;
    bool found_start_face = false;

    while (std::getline(in, line))
    {
      line = StringTrim(line);

      if (line.empty() or StartsWith(line, "//"))
        continue;

      if (line == "}")
        break;

      std::string key, value;
      if (not ReadBoundaryKeywordValue(line, key, value))
        continue;

      const auto key_l = LowerCase(key);

      if (key_l == "type")
        patch.type = value;
      else if (key_l == "nfaces")
      {
        patch.n_faces = std::stoi(value);
        found_n_faces = true;
      }
      else if (key_l == "startface")
      {
        patch.start_face = std::stoi(value);
        found_start_face = true;
      }
    }

    if (not found_n_faces or not found_start_face)
      throw std::logic_error(fname + ": Boundary patch '" + patch.name +
                             "' missing n_faces and/or start_face.");

    patches.emplace_back(std::move(patch));
  }

  return patches;
}

void
AssignBoundaryIDsFromOpenFOAMPatches(std::shared_ptr<UnpartitionedMesh> mesh,
                                     const std::vector<BoundaryPatchEntry>& patches,
                                     const std::vector<FaceLocation>& owner_face_location,
                                     size_t ninternal_faces,
                                     size_t n_faces,
                                     const std::string& fname)
{
  auto& raw_cells = mesh->GetRawCells();

  uint64_t bid = 0;
  for (const auto& patch : patches)
  {
    if (patch.n_faces < 0 or patch.start_face < 0)
      throw std::logic_error(fname + ": Negative n_faces/start_face for boundary patch '" +
                             patch.name + "'.");

    const auto start = static_cast<size_t>(patch.start_face);
    const auto end = start + static_cast<size_t>(patch.n_faces);

    if (start < ninternal_faces)
      throw std::logic_error(fname + ": Boundary patch '" + patch.name +
                             "' starts inside the internal-face range.");

    if (end > n_faces)
      throw std::logic_error(fname + ": Boundary patch '" + patch.name +
                             "' exceeds total number of faces.");

    mesh->AddBoundary(bid, patch.name);

    size_t num_assigned = 0;
    for (size_t f = start; f < end; ++f)
    {
      const auto& loc = owner_face_location.at(f);
      if (not loc.valid)
        throw std::logic_error(fname + ": Missing owner face location for boundary face " +
                               std::to_string(f) + " in patch '" + patch.name + "'.");

      auto& face = raw_cells.at(loc.cell_id)->faces.at(loc.local_face_id);

      // Adjust this to your final field name
      face.neighbor = bid;

      ++num_assigned;
    }

    log.Log() << "Assigned " << num_assigned << " faces to boundary id " << bid << "\n"
              << "\tname: " << patch.name << "\n"
              << "\ttype: " << patch.type;

    ++bid;
  }
}

std::vector<CellZoneEntry>
ReadCellZones(const std::filesystem::path& path, const std::string& fname)
{
  std::ifstream in(path);
  if (not in.is_open())
    return {}; // no cellZones is fine.
  SkipFoamHeader(in);

  // Read number of zones
  int nzones_i = 0;
  {
    std::string l;
    if (not ReadNextContentLine(in, l) or not ExtractFirstIntToken(l, nzones_i))
      throw std::logic_error(fname + ": Bad or missing zones count.");

    if (not SkipUntil(in, '('))
      throw std::logic_error(fname + ": Missing '(' starting zones list.");
  }
  const auto nzones = static_cast<size_t>(nzones_i);

  std::vector<CellZoneEntry> zones;
  zones.reserve(nzones);

  auto read_zone = [&](std::istream& in2, size_t idx) -> bool
  {
    std::string line;

    if (not ReadNextContentLine(in2, line))
      return false;
    if (line == ")")
      return false;

    CellZoneEntry z;
    z.name = line;

    if (not ReadNextContentLine(in2, line))
      return false;
    if (line != "{")
      throw std::logic_error(fname + ": Expected '{' after zone name '" + z.name + "'.");

    while (std::getline(in2, line))
    {
      line = StringTrim(line);
      if (line == "}")
        break;
      if (line.empty() or StartsWith(line, "//"))
        continue;

      if (auto semicol = line.find(';'); semicol != std::string::npos)
      {
        line.erase(semicol);
        line = StringTrim(line);
        if (line.empty())
          continue;
      }

      const std::string to_lower = LowerCase(line);
      const bool has_celllabels = to_lower.find("celllabels") != std::string::npos;
      const bool has_addressing = to_lower.find("addressing") != std::string::npos;
      const bool has_labellist = to_lower.find("labellist") != std::string::npos;

      if (has_celllabels or has_addressing or has_labellist)
      {
        int nlab = 0;

        if (not ExtractFirstIntToken(line, nlab))
        {
          std::string l2;
          while (true)
          {
            if (not ReadNextContentLine(in2, l2))
              throw std::logic_error(fname + ": EOF looking for cellLabels count in zone '" +
                                     z.name + "'.");
            if (ExtractFirstIntToken(l2, nlab) and nlab >= 0)
              break;
          }
        }

        if (not SkipUntil(in2, '('))
          throw std::logic_error(fname + ": missing '(' after count in zone '" + z.name + "'.");

        std::vector<int> labels;
        labels.reserve(static_cast<size_t>(nlab));
        if (not ReadNInts(in2, nlab, labels))
          throw std::logic_error(fname + ": Failed reading cellZones ids for zone '" + z.name +
                                 "'.");

        if (not SkipUntil(in2, ')'))
          throw std::logic_error(fname + ": missing ')' after ids in zone '" + z.name + "'.");

        z.cells = std::move(labels);
        continue;
      }
    }

    if (z.cells.empty())
      throw std::logic_error("[cellZones] Zone '" + z.name + "': no cellLabels found (empty)");

    log.Log() << "Read cellZone Id: " << idx << "\n"
              << "\tname: " << z.name << "\n"
              << "\tcell count: " << z.cells.size();

    zones.emplace_back(std::move(z));
    return true;
  };

  for (size_t i = 0; i < nzones; ++i)
  {
    if (not read_zone(in, i))
      throw std::logic_error(fname + ": Failed reading zone #" + std::to_string(i));
  }

  SkipUntil(in, ')');

  return zones;
}

} // namespace

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromOpenFOAM(const UnpartitionedMesh::Options& options)
{
  const std::string fname = "MeshIO::FromOpenFOAM";

  std::filesystem::path base = options.file_name;
  const auto poly_mesh = base / "constant" / "polyMesh";
  if (std::filesystem::is_directory(poly_mesh))
    base = poly_mesh;
  else
    throw std::logic_error(fname + ": polyMesh not found under '" + options.file_name +
                           "'. Point to <case>/constant/polyMesh or the polyMesh folder.");

  if (base.string().find("processor") != std::string::npos)
    throw std::logic_error(fname + ": Decomposed meshes (processor*/polyMesh) are not supported.");

  const auto points_p = base / "points";
  const auto faces_p = base / "faces";
  const auto owner_p = base / "owner";
  const auto neigh_p = base / "neighbour";
  const auto boundary_p = base / "boundary";

  if (not std::filesystem::exists(neigh_p))
    throw std::logic_error(fname + ": Missing 'neighbour'. Even with zero internal faces it "
                                   "should exist with count 0.");

  auto mesh = std::make_shared<UnpartitionedMesh>();
  log.Log() << "Reading OpenFOAM polyMesh from " << base.string();

  // read points

  {
    std::ifstream in(points_p);
    if (not in.is_open())
      throw std::runtime_error(fname + ": Failed to open " + points_p.string());
    SkipFoamHeader(in);

    size_t npoints = 0;
    std::vector<Vector3> verts;
    auto read_pt = [&](std::istream& in2, size_t) -> bool
    {
      double x = 0.0, y = 0.0, z = 0.0;
      if (not ReadPoint(in2, x, y, z))
        return false;
      verts.emplace_back(x, y, z);
      return true;
    };
    ReadFoamList(in, npoints, read_pt, fname, "points");

    if (verts.size() != npoints)
      throw std::logic_error(fname + ": number of points mismatch (expected " +
                             std::to_string(npoints) + ", but got " + std::to_string(verts.size()) +
                             ").");

    mesh->GetVertices() = std::move(verts);
  }

  // read faces

  std::vector<std::vector<int>> face_verts;
  {
    std::ifstream in(faces_p);
    if (not in.is_open())
      throw std::runtime_error(fname + ": Failed to open " + faces_p.string());
    SkipFoamHeader(in);

    size_t n_faces = 0;
    auto read_face = [&](std::istream& in2, size_t) -> bool
    {
      std::vector<int> fv;
      if (not ReadFace(in2, fv))
        return false;

      face_verts.emplace_back(std::move(fv));
      return true;
    };
    ReadFoamList(in, n_faces, read_face, fname, "faces");

    if (face_verts.size() != n_faces)
      throw std::logic_error(fname + ": number of faces mismatch (expected " +
                             std::to_string(n_faces) + ", but got " +
                             std::to_string(face_verts.size()) + ").");
  }

  // read owner

  std::vector<int> owner;
  {
    std::ifstream in(owner_p);
    if (not in.is_open())
      throw std::runtime_error(fname + ": Failed to open " + owner_p.string());
    SkipFoamHeader(in);

    size_t n_faces = 0;
    auto read_lab = [&](std::istream& in2, size_t) -> bool
    {
      int v = 0;
      if (not ReadFace2CellID(in2, v))
        return false;
      owner.push_back(v);
      return true;
    };
    ReadFoamList(in, n_faces, read_lab, fname, "owner");

    if (owner.size() != n_faces)
      throw std::logic_error(fname + ": owner faces count mismatch (expected " +
                             std::to_string(n_faces) + ", but got " + std::to_string(owner.size()) +
                             ").");
  }

  // read neighbour

  std::vector<int> neigh;
  {
    std::ifstream in(neigh_p);
    if (not in.is_open())
      throw std::runtime_error(fname + ": Failed to open " + neigh_p.string());
    SkipFoamHeader(in);

    size_t n_faces = 0;
    auto read_lab = [&](std::istream& in2, size_t) -> bool
    {
      int v = 0;
      if (not ReadFace2CellID(in2, v))
        return false;
      neigh.push_back(v);
      return true;
    };
    ReadFoamList(in, n_faces, read_lab, fname, "neighbour");

    if (neigh.size() != n_faces)
      throw std::logic_error(fname + ": neigh.size() != n_faces (expected " +
                             std::to_string(n_faces) + ", but got " + std::to_string(neigh.size()) +
                             ").");
  }

  // mesh count statistics

  const size_t n_faces = face_verts.size();
  const size_t ninternal_faces = neigh.size();

  if (owner.size() != n_faces)
    throw std::logic_error(fname + ": owner.size() != n_faces (expected " +
                           std::to_string(n_faces) + ", but got " + std::to_string(owner.size()) +
                           ").");

  int max_cell = -1;
  for (int v : owner)
    max_cell = std::max(max_cell, v);
  for (int v : neigh)
    max_cell = std::max(max_cell, v);
  const size_t ncells = max_cell + 1;
  if (ncells == 0)
    throw std::logic_error(fname + ": Number of cells computed cannot be zero.");

  // construct cells from the face notation in OpenFOAM

  std::vector<std::vector<int64_t>> cell_faces(ncells);

  for (size_t f = 0; f < ninternal_faces; ++f)
  {
    const int c_o = owner[f];
    const int c_n = neigh[f];

    if (c_o < 0 or std::cmp_greater_equal(c_o, ncells) or c_n < 0 or
        std::cmp_greater_equal(c_n, ncells))
      throw std::logic_error(fname + ": owner/neighbour id out of range (face " +
                             std::to_string(f) + ").");

    cell_faces[c_o].push_back(static_cast<int64_t>(f));
    cell_faces[c_n].push_back(-static_cast<int64_t>(f) - 1);
  }

  for (size_t f = ninternal_faces; f < n_faces; ++f)
  {
    const int c_o = owner[f];
    if (c_o < 0 or std::cmp_greater_equal(c_o, ncells))
      throw std::logic_error(fname + ": owner id out of range (boundary face " + std::to_string(f) +
                             ").");

    cell_faces[c_o].push_back(static_cast<int64_t>(f));
  }

  // block map from cellZones

  std::vector<unsigned int> block_map(ncells, 0);
  {
    const auto cz_path = base / "cellZones";
    const auto zones = ReadCellZones(cz_path, fname);

    if (not zones.empty())
    {
      for (size_t z_id = 0; z_id < zones.size(); ++z_id)
        for (int c_id : zones[z_id].cells)
          if (c_id >= 0 and std::cmp_less(c_id, ncells))
            block_map[c_id] = static_cast<unsigned int>(z_id);
    }
  }

  auto& raw_cells = mesh->GetRawCells();
  raw_cells.reserve(ncells);

  std::vector<FaceLocation> owner_face_location(n_faces);

  for (size_t c = 0; c < ncells; ++c)
  {
    auto cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::POLYHEDRON,
                                                                     CellType::POLYHEDRON);
    cell->block_id = block_map[c];

    for (auto code : cell_faces[c])
    {
      const int64_t f = (code >= 0) ? code : (-code - 1);
      const auto& f_v = face_verts[static_cast<size_t>(f)];

      UnpartitionedMesh::LightWeightFace lwf;
      lwf.vertex_ids.reserve(f_v.size());

      if (code >= 0)
      {
        for (int v : f_v)
          lwf.vertex_ids.push_back(static_cast<uint64_t>(v));
      }
      else
      {
        for (auto it = f_v.rbegin(); it != f_v.rend(); ++it)
          lwf.vertex_ids.push_back(static_cast<uint64_t>(*it));
      }

      const auto local_face_id = cell->faces.size();
      cell->faces.push_back(std::move(lwf));

      // Keep the owner-oriented location for each global face.
      // Boundary patch assignment will use this directly.
      if (code >= 0)
        owner_face_location.at(static_cast<size_t>(f)) = {c, local_face_id, true};
    }

    std::vector<uint64_t> v_set;
    for (const auto& f : cell->faces)
      v_set.insert(v_set.end(), f.vertex_ids.begin(), f.vertex_ids.end());

    std::sort(v_set.begin(), v_set.end());
    v_set.erase(std::unique(v_set.begin(), v_set.end()), v_set.end());
    cell->vertex_ids = std::move(v_set);

    raw_cells.push_back(std::move(cell));
  }

  mesh->SetDimension(3);
  mesh->SetType(UNSTRUCTURED);
  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  // boundary patches

  const auto patches = ReadBoundaryFile(boundary_p, fname);
  if (not patches.empty())
  {
    AssignBoundaryIDsFromOpenFOAMPatches(
      mesh, patches, owner_face_location, ninternal_faces, n_faces, fname);
  }

  log.Log() << "OpenFOAM polyMesh processed.\n"
            << "Number of nodes read: " << mesh->GetVertices().size() << "\n"
            << "Number of cells read: " << mesh->GetRawCells().size() << "\n"
            << "Number of boundary patches read: " << patches.size() << "\n";

  return mesh;
}

} // namespace opensn
