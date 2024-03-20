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
#include "utilities.hpp"
#include <queue>


//////////////////////  DFS implementation //////////////////////////////////////

size_t PlanesAlgorithmDataCleaning_DFS(OrientedTriangle T, std::vector<int> &Colors, std::vector<int> &EdgeColors, std::vector<std::vector<size_t>> &newF){

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
    std::cout << "numero de triangulos adyacentes:" << adjacent.size() << "\n";
    for(OrientedTriangle t : adjacent){
        if( Colors[t.F.getIndex()] == WHITE)
        counter += PlanesAlgorithmDataCleaning_DFS(t, Colors, EdgeColors, newF);
    }
    Colors[T.F.getIndex()] = BLACK;
    return (counter + 1);
}
/////////////////////////////////////////////////////////////////////////

//////////////////////////  BFS implementation /////////////////////////

size_t PlanesAlgorithmDataCleaning_BFS(OrientedTriangle T0, std::vector<int> &Colors, std::vector<int> &EdgeColors, std::vector<std::vector<size_t>> &newF){

    std::queue<OrientedTriangle> Queue;
    Queue.push(T0);
    size_t counter = 0;
    std::vector<OrientedTriangle> adjacent; // list of adjacent triangles of T that shall be explored 

    while(!Queue.empty()){
        counter++;
        OrientedTriangle T = Queue.front();
        Queue.pop();
        std::cout << "Triangle: " << T.F.getIndex() << "\n";

        Vector3 N1 = normal_Vec(T); // normal vector of triangle T
        std::vector<size_t> triangle = {T.v1, T.v2, T.v3}; 
        newF.push_back(triangle); // add to the list of new triangles
        adjacent.clear();

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

        std::cout << "numero de triangulos adyacentes:" << adjacent.size() << "\n";
        for(OrientedTriangle t : adjacent){
            if( Colors[t.F.getIndex()] == WHITE){
                Queue.push(t);
                Colors[t.F.getIndex()] = GRAY;
            }
        }
        Colors[T.F.getIndex()] = BLACK;
    }
    return counter;
}
/////////////////////////////////////////////////////////////////////////
