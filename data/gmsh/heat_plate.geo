SetFactory("OpenCASCADE");

// create circular domain
Disk(1) = {0, 0, 0, 2.2, 2.2};

// set mesh size
MeshSize{:} = 0.2;

// create and save mesh 
Mesh 2;
Save "heat_plate.msh";
