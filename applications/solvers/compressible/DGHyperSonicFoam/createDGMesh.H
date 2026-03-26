// dgGeomMesh is the bridge between fvMesh and DG operators. It owns the
// DG cells/faces, precomputed quadrature data, basis information, and cell
// metrics shared by all DG fields and DG operators. A useful mental model is
// that fvMesh gives the topological mesh, while dgGeomMesh adds the machinery
// needed to evaluate/interpolate/integrate high-order DG fields on that mesh.
dgGeomMesh dgMesh(mesh);

Info<< "Create DG Geometric Mesh with polynomial order "
    << dgMesh.pOrder() << "\n" << endl;

// dgMesh.cells()[0]->printDebugInfo();
