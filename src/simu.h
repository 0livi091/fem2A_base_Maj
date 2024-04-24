#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            //std::cout << "TO BE IMPLEMENTED !!!" << std::endl;
            
            Mesh mesh;
            mesh.load(mesh_filename);
            Quadrature quad = Quadrature::get_quadrature(2);
            SparseMatrix K_glob(mesh.nb_vertices());
            ShapeFunctions fonctions(2,1);
            
            std::vector<double> F_glob(mesh.nb_vertices(), 0.);
            
            mesh.set_attribute(unit_fct, 1, true);
            
            std::vector <bool> attribute_is_dirichlet(2, false);
            attribute_is_dirichlet[1]=true;
            
            
            std::vector<double> values(mesh.nb_vertices());
            for(int e=0; e<mesh.nb_vertices(); e++){
            	values[e] = xy_fct(mesh.get_vertex(e));
            }
            
            for(int t=0 ;t<mesh.nb_triangles() ; t++){
    
            	
            	ElementMapping elt(mesh, false, t);
            	DenseMatrix Ke;
            	Ke.set_size(3,3);
            	assemble_elementary_matrix(elt, fonctions, quad, unit_fct, Ke);
     
       	    	local_to_global_matrix(mesh, t, Ke, K_glob);
       	    	
       	    }
       	    
       	    apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet,values,K_glob,F_glob);
       	    std::vector<double> u(mesh.nb_vertices());
       	    solve(K_glob, F_glob, u);
       	    std::string export_name = "pure_dirichlet_square_fine";
       	    mesh.save(export_name+".mesh");
       	    save_solution(u, export_name+".bb");
        }


// avec terme source :
// values vecteur vaut 0
// terme source non nul avec vecteur F 


	void dirichlet_terme_source_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            //std::cout << "TO BE IMPLEMENTED !!!" << std::endl;
            
            Mesh mesh;
            mesh.load(mesh_filename);
            Quadrature quad = Quadrature::get_quadrature(2);
            SparseMatrix K_glob(mesh.nb_vertices());
            ShapeFunctions fonctions(2,1);
            
            std::vector< double > Fe(fonctions.nb_functions(), 0.0);
            std::vector<double> F_glob(mesh.nb_vertices(), 0.);
            
            mesh.set_attribute(unit_fct, 1, true);
            
            std::vector <bool> attribute_is_dirichlet(2, false);
            attribute_is_dirichlet[1]=true;
            
            
            std::vector<double> values(mesh.nb_vertices(), 0.0);
            
            for(int t=0 ;t<mesh.nb_triangles() ; t++){
    
            	
            	ElementMapping elt(mesh, false, t);
            	DenseMatrix Ke;
            	Ke.set_size(3,3);
            	assemble_elementary_matrix(elt, fonctions, quad, unit_fct, Ke);
     
       	    	local_to_global_matrix(mesh, t, Ke, K_glob);
       	    	
            	assemble_elementary_vector(elt, fonctions, quad, unit_fct, Fe);
           
            	local_to_global_vector(mesh, false, t, Fe, F_glob);
       	    	
       	    }
       	    
       	    apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet,values,K_glob,F_glob);
       	    std::vector<double> u(mesh.nb_vertices());
       	    
       	    solve(K_glob, F_glob, u);
       	    
       	    std::string export_name = "dirichlet_sourceterm_square";
       	    mesh.save(export_name+".mesh");
       	    save_solution(u, export_name+".bb");
        }
     }
}
