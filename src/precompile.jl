using PrecompileTools: @compile_workload

@compile_workload begin
    m = MMJMesh.Gmsh.Mesh("""
        SetFactory("OpenCASCADE");
        Rectangle(1) = {0, 0, 0, 20, 10, 0};
        Rectangle(2) = {10, -5, 0, 5, 20, 0};
        Disk(3) = {20, 5, 0, 15, 5};
        BooleanUnion{
            Surface{1}; 
            Delete;
        }{
            Surface{2}; 
            Surface{3}; 
            Delete;
        }
        MeshSize{:} = 0.5;
        Mesh 2;
    """
    )
    MMJMesh.Plots.mplot(m, edgesvisible=true) |> MMJMesh.Plots.mconf()
    MMJMesh.Plots.mplot(m, rand(nedges(m)), edgesvisible=true) |> MMJMesh.Plots.mconf()
    MMJMesh.Plots.mplot(m, rand(nfaces(m)), edgesvisible=true) |> MMJMesh.Plots.mconf()
    MMJMesh.Plots.mplot(m, rand(2 * nedges(m)), edgesvisible=true) |> MMJMesh.Plots.mconf()
end
