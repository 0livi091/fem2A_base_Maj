#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
       bool test_quadrature(int order)
       {
            Quadrature quad = Quadrature::get_quadrature(order);
            std :: cout << quad.nb_points() << std ::endl;
            double sum = 0;
              for (int i=0; i< quad.nb_points();++i){
              std::cout << quad.point(i).x << " " // point a été défini comme une structure vertex (dans mesh.h)
              << quad.point(i).y << std::endl;
              std::cout << quad.weight(i) << std::endl;
              sum = sum + quad.weight(i);
              }
            std::cout << sum << std::endl;
            return true;
        }
        
       bool test_elementmapping(){
       
            std :: cout << "Coordonnées des éléments : " << std :: endl;
       
            Mesh mesh;
            mesh.load("data/square.mesh");
       
            ElementMapping e4(mesh, true, 4);
            ElementMapping e42(mesh, false, 4);
            
            return true;
            
       }
       
       bool test_elementtransform(){
       
            std :: cout << "Mapping : \n" ;
            
            Mesh mesh;
            mesh.load("data/square.mesh");
            
            //ElementMapping e4(mesh, true, 4);
            ElementMapping e42(mesh, false, 4);
            
            vertex vr ;
            vr.x = 0.2;
            vr.y = 0.4;
            
            //vertex vr_new = e4.transform(vr);
            vertex vr_new2 = e42.transform(vr);
            
            //std:: cout << "x = " << vr_new.x << " ; y = " << vr_new.y << std :: endl;
            std:: cout << "x = " << vr_new2.x << " ; y = " << vr_new2.y << std :: endl;
            
            return true;
       }
       

    }
}
