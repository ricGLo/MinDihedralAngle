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

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

////////////////////////////////////////     Planes Algorithm     ///////////////////////////////////////////

///////////////  Utilities
std::vector<int> Colors; // Face Colors
std::vector<int> EdgeColors;  // Edge Colors
std::vector<std::vector<size_t>> newF;  // List of newFaces 
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

// This function computes the normal vector of a OrientedTriangle 
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

//////////////////////  Implementation    

size_t PlanesAlgorithmDataCleaning_DFS(OrientedTriangle T){

  Vector3 N1 = normal_Vec(T); // normal vector of triangle T
  std::vector<size_t> triangle = {T.v1, T.v2, T.v3}; 
  newF.push_back(triangle); // add to the list of new triangles 
  std::vector<OrientedTriangle> adjacent; // list of adjacent triangles of T that shall be explored
  size_t counter = 0;


  for(Edge e: T.F.adjacentEdges()){ // Fix manifoldness of each edge

    if( EdgeColors[e.getIndex()] == WHITE ){  // If edge e hasn't been explored

    OrientedTriangle T_prime = ReOrder(T.F, e, T.v1, T.v2, T.v3); // T_prime is the same T and has same orientation


      size_t bestP;
      Face bestF;
      double alpha_min = 100.0;
      for(Face f: e.adjacentFaces()){ // Find the face f adjacent to edge e whose dihedral angle with T is minimized 

        if(f.getIndex() != T.F.getIndex()){  // If f is different from T 

          size_t p = get_third(T_prime.v1, T_prime.v2, f); // get the vertex of f which is different from v1 and v2
              // get normal vector of face f with the same relative orientation as T
          Vector3 N2 = normal_Vec(OrientedTriangle(f, T_prime.v2, T_prime.v1, p)); 
          double alpha = angle(N1, N2);
          if(alpha < alpha_min){
            alpha_min = alpha;
            bestP = p;
            bestF = f;
          }    
        }   
      }
      if(!e.isBoundary()) // if e is not boundary it is part of at least two triangles 
        adjacent.push_back(OrientedTriangle(bestF, T_prime.v2, T_prime.v1, bestP));  

      EdgeColors[e.getIndex()] = BLACK; // mark as explored 
    }
    
  }

  Colors[T.F.getIndex()] = GRAY;  // mark T as explored 
  for(OrientedTriangle t : adjacent){
    if( Colors[t.F.getIndex()] == WHITE)
      counter += PlanesAlgorithmDataCleaning_DFS(t);
  }
  Colors[T.F.getIndex()] = BLACK;
  return (counter + 1);
}
//////////////////////////////////////////////////


int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();

  // Load mesh
  std::tie(mesh, geometry) = readSurfaceMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(args::get(inputFilename)),
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));


  //////////////// CODE HERE  ////////////////////

  // Assign memory 
  EdgeColors.assign(mesh->nEdges(), WHITE);
  Colors.assign(mesh->nFaces(), WHITE);


  Face F = mesh->face(0);

  size_t v1 = 0;
  size_t v2 = 2;
  size_t v3 = 1;

  OrientedTriangle T(F, 0, 2, 1);

  size_t C = PlanesAlgorithmDataCleaning_DFS(T);

  std::cout << "Triangles found = " << C << "\n\n";

  for(auto t : newF){
    for(size_t i : t)
      std::cout << i+1 << "\t";
    std::cout<< "\n";
  }
  ////////////////////////////////////////

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
