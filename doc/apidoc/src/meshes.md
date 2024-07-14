# Meshes

```@meta
CurrentModule = MMJMesh.Meshes
```

The `MMJMesh.Meshes` module is the main interface to the functionality of `MMJMesh`. Its main datastructure is `Mesh`

```@docs
Mesh
```

Example:

```@example
using MMJMesh.Meshes
Mesh(
    [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9], 
    [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]], 
    2
)
```
