// SPDX-FileCopyrightText: 2024 The OpenSn Authors
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

// system includes
#include <algorithm>
#include <fstream>

namespace opensn
{
namespace
{
// Utility functions

// trim whitespace from the string
void
TrimInPlace(std::string& str)
{
  auto is_space = [](unsigned char ch) { return std::isspace(ch) != 0; };

  auto first = std::find_if_not(
    str.begin(), str.end(), [&](char ch) { return is_space(static_cast<unsigned char>(ch)); });
  auto last = std::find_if_not(str.rbegin(),
                               str.rend(),
                               [&](char ch) { return is_space(static_cast<unsigned char>(ch)); })
                .base();

  if (first >= last)
  {
    str.clear();
    return;
  }

  str.erase(last, str.end());
  str.erase(str.begin(), first);
}

// verify file section identifiers
inline bool
StartsWithKey(std::string_view str, std::string_view key)
{
  return str.size() >= key.size() && std::equal(key.begin(), key.end(), str.begin());
}

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
      while (s && !(prev == '*' && (cur = s.get()) == '/'))
        prev = cur;
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

  if (!(in >> x >> y >> z))
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
  if (!(in >> n) || n < 0)
    return false;

  // find opening '('
  if (!ExpectChar(in, '('))
    return false;

  // read n vertex labels
  out.clear();
  out.reserve(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
  {
    SkipWSAndComments(in);
    int v;
    if (!(in >> v))
      return false;
    out.push_back(v);
  }

  // find closing ')'
  if (!ExpectChar(in, ')'))
    return false;

  return true;
}

// 'owner'/'neighbour' files
inline bool
ReadLabel(std::istream& in, int& v)
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
    std::string tmp = line;
    TrimInPlace(tmp);
    if (tmp.empty() || StartsWithKey(tmp, "//"))
      continue;
    if (!tmp.empty() && (std::isdigit(static_cast<unsigned char>(tmp[0])) || tmp[0] == '-'))
    {
      // rewind to start of this line
      in.seekg(-static_cast<std::streamoff>(line.size()) - 1, std::ios_base::cur);
      return;
    }
  }
  in.clear();
}

// skip chars until either a ' ( ' or ' ) ' is encountered
inline bool
SkipUntil(std::istream& in, char target)
{
  int ch;
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
    int v;
    if (!(in >> v))
      return false;
    out.emplace_back(v);
  }
  return true;
}

struct CellZoneEntry
{
  std::string name;
  std::vector<int> cells;
};

// generalized function to parse/select OpenFOAM Lists (data structure in file)
// Supports both paren lists:  N ( item0 item1 ... )
// and uniform brace lists:    N { value }   ==> repeats 'value' N times.
template <class ItemReader>
void
ReadFoamList(
  std::istream& in, int& count_out, ItemReader item_fn, const std::string& fname, const char* name)
{

  // read count
  SkipWSAndComments(in);
  if (!(in >> count_out) || count_out < 0)
    throw std::logic_error(fname + ": Bad " + std::string(name) + " count.");

  // detect opener: '(' ( list ) or '{' ( uniform repeat )
  SkipWSAndComments(in);
  int opener = in.peek();
  if (opener != '(' && opener != '{')
    throw std::logic_error(fname + ": Expected '(' or '{' after " + std::string(name) + " count.");
  in.get(); // consume opener

  const bool brace_uniform = (opener == '{');

  if (!brace_uniform)
  {
    // Parenthesis list: read count_out items normally
    for (int i = 0; i < count_out; ++i)
      if (!item_fn(in, i))
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
        if (!in)
          break;
        if (ch == ')')
        {
          closed = true;
          break;
        }
      }
      if (!closed)
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
    for (int i = 0; i < count_out; ++i)
    {
      in.clear();
      in.seekg(value_pos);
      if (!item_fn(in, i))
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
        if (!in)
          break;
        if (ch == '}')
        {
          closed = true;
          break;
        }
      }
      if (!closed)
        throw std::logic_error(fname + ": Missing '}' for " + std::string(name));
    }
  }
}

std::vector<CellZoneEntry>
ReadCellZones(const std::filesystem::path& path, const std::string& fname)
{
  std::ifstream in(path);
  if (!in.is_open())
    return {}; // no cellZones is fine.
  SkipFoamHeader(in);

  int nzones = 0;
  std::vector<CellZoneEntry> zones;

  auto read_zone = [&](std::istream& in2, int /*idx*/) -> bool
  {
    std::string line;

    // read zone name
    if (!std::getline(in2, line))
      return false;
    TrimInPlace(line);
    if (line.empty())
      return false;

    CellZoneEntry z;
    z.name = line;

    // validity check
    if (!std::getline(in2, line))
      return false;
    TrimInPlace(line);
    if (line != "{")
      return false;

    // read buffer
    while (std::getline(in2, line))
    {
      TrimInPlace(line);
      if (line == "}")
        break;
      if (line.empty() || StartsWithKey(line, "//"))
        continue;

      auto semicol = line.find(';');
      if (semicol != std::string::npos)
        line.erase(semicol);
      TrimInPlace(line);

      if (StartsWithKey(line, "cellLabels") || StartsWithKey(line, "addressing") ||
          StartsWithKey(line, "labelList"))
      {
        // next non-empty content should be the count, then a parenthesis block with IDs.

        // read the count entry
        int nlab = 0;
        {
          std::streampos pos_before;
          std::string l2;
          while (true)
          {
            pos_before = in2.tellg();
            if (!std::getline(in2, l2))
              return false;
            TrimInPlace(l2);
            if (l2.empty() || StartsWithKey(l2, "//"))
              continue;
            // Put back so we can read the integer via operator>>
            in2.clear();
            in2.seekg(pos_before);
            break;
          }
          if (!(in2 >> nlab) || nlab < 0)
            throw std::logic_error(fname + ": Bad cellZones count for cellLabels/addressing.");
        }

        // start parse at '('
        if (!SkipUntil(in2, '('))
          throw std::logic_error(fname + ": cellZones: missing '(' after cellLabels count.");

        // read through the cellZones
        std::vector<int> labels;
        if (!ReadNInts(in2, nlab, labels))
          throw std::logic_error(fname + ": Failed reading cellZones ids.");

        // parse until section until ')'
        if (!SkipUntil(in2, ')'))
          throw std::logic_error(fname + ": cellZones: missing ')' after ids.");

        z.cells = std::move(labels);
      }
    }

    zones.emplace_back(z);
    return true;
  };

  ReadFoamList(in, nzones, read_zone, fname, "cellZones");
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

  if (!std::filesystem::exists(neigh_p))
    throw std::logic_error(fname + ": Missing 'neighbour'. Even with zero internal faces it "
                                   "should exist with count 0.");

  auto mesh = std::make_shared<UnpartitionedMesh>();
  log.Log() << "Reading OpenFOAM polyMesh from " << base.string();

  // read " points "
  {
    std::ifstream in(points_p);
    if (!in.is_open())
      throw std::runtime_error(fname + ": Failed to open " + points_p.string());
    SkipFoamHeader(in);

    int npoints = 0;
    std::vector<Vector3> verts;
    auto read_pt = [&](std::istream& in2, int) -> bool
    {
      double x, y, z;
      if (!ReadPoint(in2, x, y, z))
        return false;
      verts.emplace_back(x, y, z);
      return true;
    };
    ReadFoamList(in, npoints, read_pt, fname, "points");
    if (verts.size() != npoints)
    {
      throw std::logic_error(fname + ": number of points mismatch (expected " +
                             std::to_string(npoints) + ", but got " + std::to_string(verts.size()) +
                             ").");
    }
    mesh->GetVertices() = std::move(verts);
  }

  // read " faces "
  std::vector<std::vector<int>> face_verts;
  {
    std::ifstream in(faces_p);
    if (!in.is_open())
      throw std::runtime_error(fname + ": Failed to open " + faces_p.string());
    SkipFoamHeader(in);

    int nfaces = 0;
    auto read_face = [&](std::istream& in2, int) -> bool
    {
      std::vector<int> fv;
      if (!ReadFace(in2, fv))
      {
        return false;
      }
      else
      {
        face_verts.emplace_back(std::move(fv));
        return true;
      }
    };
    ReadFoamList(in, nfaces, read_face, fname, "faces");
    if (face_verts.size() != nfaces)
    {
      throw std::logic_error(fname + " : number of faces mismatch (expected " +
                             std::to_string(nfaces) + ", but got " +
                             std::to_string(face_verts.size()) + ").");
    }
  }

  // read " owner "
  std::vector<int> owner;
  {
    std::ifstream in(owner_p);
    if (!in.is_open())
      throw std::runtime_error(fname + ": Failed to open " + owner_p.string());
    SkipFoamHeader(in);

    int nfaces = 0;
    auto read_lab = [&](std::istream& in2, int) -> bool
    {
      int v;
      if (!ReadLabel(in2, v))
        return false;
      owner.push_back(v);
      return true;
    };
    ReadFoamList(in, nfaces, read_lab, fname, "owner");
    if (owner.size() != nfaces)
    {
      throw std::logic_error(fname + ": owner faces count mismatch (expected " +
                             std::to_string(nfaces) + ", but got " + std::to_string(owner.size()) +
                             ").");
    }
  }

  // read " neighbour " ( must exist, but can have a zero entry )
  std::vector<int> neigh;
  {
    std::ifstream in(neigh_p);
    if (!in.is_open())
      throw std::runtime_error(fname + ": Failed to open " + neigh_p.string());
    SkipFoamHeader(in);

    int nfaces = 0;
    auto read_lab = [&](std::istream& in2, int) -> bool
    {
      int v;
      if (!ReadLabel(in2, v))
        return false;
      neigh.push_back(v);
      return true;
    };
    ReadFoamList(in, nfaces, read_lab, fname, "neighbour");
    if (neigh.size() != nfaces)
    {
      throw std::logic_error(fname + ": neigh.size() != nfaces ( expected " +
                             std::to_string(nfaces) + ", but got " + std::to_string(neigh.size()) +
                             ").");
    }
  }

  // construct cells
  const int nfaces = static_cast<int>(face_verts.size());
  const int ninternal_faces = static_cast<int>(neigh.size());
  if (owner.size() != nfaces)
  {
    throw std::logic_error(fname + ": owner.size() != nfaces ( expected " + std::to_string(nfaces) +
                           ", but got " + std::to_string(owner.size()) + ").");
  }

  int max_cell = -1;
  for (int v : owner)
    max_cell = std::max(max_cell, v);
  for (int v : neigh)
    max_cell = std::max(max_cell, v);
  const int ncells = max_cell + 1;
  if (ncells <= 0)
    throw std::logic_error(fname + ": Non-positive number of cells computed.");

  auto& raw_cells = mesh->GetRawCells();
  raw_cells.reserve(ncells);

  // for each cell, collect faces according to OpenFOAM orientation
  std::vector<std::vector<int>> cell_faces(ncells);
  for (int f = 0; f < ninternal_faces; ++f)
  {
    const int c_o = owner[f];
    const int c_n = neigh[f];
    if (c_o < 0 || c_o >= ncells || c_n < 0 || c_n >= ncells)
      throw std::logic_error(fname + ": owner/neighbour id out of range (face " +
                             std::to_string(f) + ").");
    cell_faces[c_o].push_back(+f);     // owner side: outward as stored
    cell_faces[c_n].push_back(-f - 1); // neighbour side: reversed for outward
  }
  for (int f = ninternal_faces; f < nfaces; ++f)
  {
    const int c_o = owner[f];
    if (c_o < 0 || c_o >= ncells)
      throw std::logic_error(fname + ": owner id out of range (boundary face " + std::to_string(f) +
                             ").");
    cell_faces[c_o].push_back(+f); // boundary face outward w.r.t owner
  }

  // volume cells (POLYHEDRON)
  for (int c = 0; c < ncells; ++c)
  {
    auto cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::POLYHEDRON,
                                                                     CellType::POLYHEDRON);
    cell->block_id = 0;

    for (int code : cell_faces[c])
    {
      const int f = (code >= 0) ? code : (-code - 1);
      const auto& f_v = face_verts[f];

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
      cell->faces.push_back(std::move(lwf));
    }

    // per cell vertex ids
    std::vector<uint64_t> v_set;
    for (const auto& f : cell->faces)
      v_set.insert(v_set.end(), f.vertex_ids.begin(), f.vertex_ids.end());
    std::sort(v_set.begin(), v_set.end());
    v_set.erase(std::unique(v_set.begin(), v_set.end()), v_set.end());
    cell->vertex_ids = std::move(v_set);

    raw_cells.push_back(std::move(cell));
  }

  // assign block ids
  const auto cz_path = base / "cellZones";
  const auto zones = ReadCellZones(cz_path, fname);

  if (!zones.empty())
  {
    for (int z_id = 0; z_id < (int)zones.size(); ++z_id)
    {
      for (int c_id : zones[z_id].cells)
      {
        if (c_id >= 0 && c_id < (int)raw_cells.size())
          raw_cells[c_id]->block_id = z_id;
      }
    }
  }

  mesh->SetDimension(3);
  mesh->SetType(UNSTRUCTURED);
  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  log.Log() << "  OpenFOAM polyMesh processed.\n"
            << "  Number of nodes read     : " << mesh->GetVertices().size() << "\n"
            << "  Number of cells read  : " << mesh->GetRawCells().size() << "\n";

  return mesh;
}

} // namespace opensn
