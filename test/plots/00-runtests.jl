module PlotsTest

using Test
import CairoMakie

using MMJMesh
using MMJMesh.Plots
using MMJMesh.Meshes
using MMJMesh.Utilities

@testset "Helperfunctions" include("helperfunctionstests.jl")
@testset "Sample1D" include("sample1dtests.jl")
@testset "CurveApproximation" include("curveapproximationtests.jl")
@testset "Reference" include("referencetests.jl")

end
