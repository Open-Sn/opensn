// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

SetFactory("OpenCASCADE");
Box(1) = {-30, -30, -30, 120, 160, 120};
Mesh.MeshSizeMin = 20;
Mesh.MeshSizeMax = 20;
Mesh 2;
Save "kobayashi_box.stl";
