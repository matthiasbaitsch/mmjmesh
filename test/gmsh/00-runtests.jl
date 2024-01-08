# module GmshTests

using Test
using MMJMesh
using MMJMesh.Gmsh
using MMJMesh.Meshes

# @testset "Gmsh" begin

    # Path to input file
    fp(filename) = dirname(@__FILE__) * "/../../data/gmsh/" * filename
    
    # Test that files are processed without error
    m = MMJMesh.Gmsh.readmesh(fp("advanced.msh"))
    @test m.nodeBlocks.nnodes == 21

    # m = MMJMesh.Gmsh.Mesh(fp("advanced.msh"))
    # @test nnodes(m) == 21

    # XXX This one fails
    # m = MMJMesh.Gmsh.readmesh(fp("simple.msh"))    
# end

# end