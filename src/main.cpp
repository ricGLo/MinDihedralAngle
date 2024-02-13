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

///  Utilities
std::vector<int> Colors;
std::vector<int> EdgeColors;
std::vector<std::vector<size_t>> newF;
                    //MACROS
#define BLACK 2
#define WHITE 0
#define GRAY 1

struct OrientedTriangle{
  Face F;
  // Vertices 
  size_t v1, v2, v3;

  OrientedTriangle(Face F_, size_t v1_, size_t v2_, size_t v3_){
    F = F_; v1 = v1_; v2 = v2_, v3 = v3_;
  }
};

Vector3 normal_Vec(OrientedTriangle T){
  Vector3 V1 = geometry->vertexPositions[mesh->vertex(T.v1)];
  Vector3 V2 = geometry->vertexPositions[mesh->vertex(T.v2)];
  Vector3 V3 = geometry->vertexPositions[mesh->vertex(T.v3)];

  Vector3 N = cross(V2-V1, V3-V1);
  return unit(N);
}

OrientedTriangle ReOrder(Face f, Edge e, size_t v1, size_t v2, size_t v3){
  std::array<Vertex, 2> verts = e.adjacentVertices();

  if( (v1 == verts[0].getIndex() and v2 == verts[1].getIndex()) or (v2 == verts[0].getIndex() and v1 == verts[1].getIndex()) ) 
    return OrientedTriangle(f, v1, v2, v3);
  else if( (v1 == verts[0].getIndex() and v3 == verts[1].getIndex()) or (v3 == verts[0].getIndex() and v1 == verts[1].getIndex()) )
    return OrientedTriangle(f, v3, v1, v2);
  else 
    return OrientedTriangle(f, v2, v3, v1);
}

size_t get_third(size_t v1, size_t v2, Face f){
  for( Vertex vert: f.adjacentVertices()){
    size_t idx = vert.getIndex();
    if (idx != v1 and idx != v2)
      return idx;
  }
}

//////////////////////  Implementation    

size_t PlanesAlgorithmDataCleaning(OrientedTriangle T, size_t counter){
  Vector3 N1 = normal_Vec(T);
  std::vector<size_t> triangle = {T.v1, T.v2, T.v3};
  newF.push_back(triangle);
  std::vector<OrientedTriangle> adjacent;

  for(Edge e: T.F.adjacentEdges()){
    OrientedTriangle T_prime = ReOrder(T.F, e, T.v1, T.v2, T.v3);

    if( EdgeColors[e.getIndex()] == WHITE ){

      size_t bestP;
      Face bestF;
      double alpha_min = 100.0;
      for(Face f: e.adjacentFaces()){

        if(f.getIndex() != T.F.getIndex()){  // If F is different from T 

          size_t p = get_third(T_prime.v1, T_prime.v2, f); // get the vertex of f which is different from v1 and v2
          Vector3 N2 = normal_Vec(OrientedTriangle(f, T_prime.v2, T_prime.v1, p));
          double alpha = angle(N1, N2);
          if(alpha < alpha_min){
            alpha_min = alpha;
            bestP = p;
            bestF = f;
          }    
        }   
      }
      if(!e.isBoundary())
        adjacent.push_back(OrientedTriangle(bestF, T_prime.v2, T_prime.v1, bestP));  
    }
    EdgeColors[e.getIndex()] = BLACK;
  }

  Colors[T.F.getIndex()] = GRAY;
  for(OrientedTriangle t : adjacent){
    if( Colors[t.F.getIndex()] == WHITE)
      counter += PlanesAlgorithmDataCleaning(t, counter);
  }
  Colors[T.F.getIndex()] = BLACK;
  return (counter + 1);
}


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


  //////////////// TODO ////////////////////

  Face F = mesh->face(0);

  size_t v1 = 0;
  size_t v2 = 2;
  size_t v3 = 1;

  OrientedTriangle T(F, 0, 2, 1);

  size_t C = PlanesAlgorithmDataCleaning(T, 0);

  // std::cout << "Triangles found = " << C << "\n\n";

  for(auto t : newF){
    for(size_t i : t)
      std::cout << i << "\t";
    std::cout<< "\n";
  }


  ////////////////////////////////////////

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
