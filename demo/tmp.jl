using GeometryBasics, LinearAlgebra, GLMakie, FileIO



#mesh(m)


brain = load(assetpath("brain.stl"))

mesh(brain)