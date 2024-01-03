import CairoMakie as cm

using MMJMesh
using MMJMesh.Plots
using MMJMesh.Meshes
using MMJMesh.Utilities
using MMJMesh.Topologies
using MMJMesh.Mathematics

cm.set_theme!(cm.theme_minimal())
cm.update_theme!(colormap=:jet)