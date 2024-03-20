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
#include "MinDihedral.hpp"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;


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

  std::vector<std::vector<size_t>> newF;  // List of newFaces 
  std::vector<int> EdgeColors;
  std::vector<int> Colors;
  // Assign memory 
  EdgeColors.assign(mesh->nEdges(), WHITE);
  Colors.assign(mesh->nFaces(), WHITE);

  Face F = mesh->face(0);

  size_t v1 = 0;
  size_t v2 = 2;
  size_t v3 = 1;

  OrientedTriangle T(F, 0, 2, 1);

  size_t C = PlanesAlgorithmDataCleaning_BFS(T, Colors, EdgeColors, newF);

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
