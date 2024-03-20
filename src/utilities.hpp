#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/direction_fields.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/surface_mesh.h" 
#include "args/args.hxx"
#include "imgui.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
ManifoldSurfaceMesh FinalMesh;


///////////////  Utilities

//Colors
#define BLACK 2 // visited
#define WHITE 0 // un-explored
#define GRAY 1 // explored

// Structure of an oriented face. This allows to orient a face oppositely as in the original mesh
struct OrientedTriangle{
  // face
  Face F;
  // Vertices 
  size_t v1, v2, v3;

  //constructor
  OrientedTriangle(Face F_, size_t v1_, size_t v2_, size_t v3_){
    F = F_; v1 = v1_; v2 = v2_, v3 = v3_;
  }
};

// This function computes the normal vector of an OrientedTriangle 
Vector3 normal_Vec(OrientedTriangle T){
  Vector3 V1 = geometry->vertexPositions[mesh->vertex(T.v1)];
  Vector3 V2 = geometry->vertexPositions[mesh->vertex(T.v2)];
  Vector3 V3 = geometry->vertexPositions[mesh->vertex(T.v3)];

  Vector3 N = cross(V2-V1, V3-V1);
  return unit(N);
}


// INPUT: face f, edge e, vertices v1, v2, v3
// OUTPUT: an OrientedTriangle corresponding to face f with orientation u1, u2, u3 where [u1, u2] is the edge e
//          and such that u1, u2, u3 is a even permutation of v1, v2, v3
OrientedTriangle ReOrder(Face f, Edge e, size_t v1, size_t v2, size_t v3){
  std::array<Vertex, 2> verts = e.adjacentVertices();

  if( (v1 == verts[0].getIndex() and v2 == verts[1].getIndex()) or (v2 == verts[0].getIndex() and v1 == verts[1].getIndex()) ) 
    return OrientedTriangle(f, v1, v2, v3);
  else if( (v1 == verts[0].getIndex() and v3 == verts[1].getIndex()) or (v3 == verts[0].getIndex() and v1 == verts[1].getIndex()) )
    return OrientedTriangle(f, v3, v1, v2);
  else 
    return OrientedTriangle(f, v2, v3, v1);
}

// Return the vertex of f different from v1 and v2
size_t get_third(size_t v1, size_t v2, Face f){
  for( Vertex vert: f.adjacentVertices()){
    size_t idx = vert.getIndex();
    if (idx != v1 and idx != v2)
      return idx;
  }
  return INVALID_IND;
}