#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border )
    {
    // border true si on est sur un bord = un segment avec 2 coordonnees = 2vertices
    // border false = triangle, 3 coordonnees, 3 vertices
    // v est local, i est global (numero du triangle)
    	if (border){
    		for (int v=0; v<2 ; ++v) vertices_.push_back(M.get_edge_vertex(i,v));
    		//for (int v=0; v<2 ; ++v) std::cout << vertices_[v].x << " "<< vertices_[v].y;
    	}
    	else {
    		for (int v=0; v<3; ++v) vertices_.push_back(M.get_triangle_vertex(i,v));
    		//for (int v=0; v<3; ++v) std::cout << vertices_[v].x << " "<< vertices_[v].y << std::endl;
    	} 
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
    
    	double xi = x_r.x;
    	double eta = x_r.y;
    	double sum = 0;
    	vertex r ;
    	
    	if (border_){
    	
    		r.x = vertices_[0].x*(1-xi) + vertices_[1].x*xi;
    		r.y = vertices_[0].y*(1-xi) + vertices_[1].y*xi;
    	}
    	else {
    	
    		r.x = vertices_[0].x*(1-xi-eta) + vertices_[1].x*xi + vertices_[2].x*eta;
    		r.y = vertices_[0].y*(1-xi-eta) + vertices_[1].y*xi + vertices_[2].y*eta;
    	}
    	
    	
        //std::cout << "[ElementMapping] transform reference to world space" << '\n';

        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] compute jacobian matrix" << '\n';
        DenseMatrix J ;
        
        double xi = x_r.x;
        double eta = x_r.y;
        
        if(border_){
        	J.set_size(2,1);
        	J.set(0,0,-vertices_[0].x + vertices_[1].x);
        	J.set(1,0,-vertices_[0].y + vertices_[1].y);
        }
        
        else {
        	J.set_size(2,2);
        	J.set(0,0,-vertices_[0].x + vertices_[1].x);
        	J.set(0,1,-vertices_[0].x + vertices_[2].x);
        	J.set(1,0,-vertices_[0].y + vertices_[1].y);
        	J.set(1,1,-vertices_[0].y + vertices_[2].y);
        }
        
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const 
    {
        //std::cout << "[ElementMapping] compute jacobian determinant" << '\n';
        DenseMatrix jacob = jacobian_matrix(x_r);
        double det = 0;
        if(border_){
        	double det = std::sqrt(jacob.get(0,0)+jacob.get(0,0) + jacob.get(1,0)+jacob.get(1,0));
        }
        else{
        	double det = jacob.det_2x2();
        }
        return det ;
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        //std::cout << "[ShapeFunctions] constructor in dimension " << dim << '\n';
        bool SF_construct = true;
        if( dim!= 1 && dim!= 2){
        	std::cout << "ShapeFunctions are only implemented in 1D or 2D. \n";
        	SF_construct = false;
        }
        if (order != 1){
        	std::cout << "Only order-1 ShapeFunctions are implemented. \n";
        	SF_construct = false;
        }
        assert(SF_construct);

    }

    int ShapeFunctions::nb_functions() const
    {
        //std::cout << "[ShapeFunctions] number of functions" << '\n';
        if (dim_ == 1) return 2 ;
        return 3;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        //std::cout << "[ShapeFunctions] evaluate shape function " << i << '\n';
        if (dim_ == 1){
        	switch (i){
        		case (0) : 
        			return 1 - x_r.x;
        		case (1) : 
        			return x_r.x;
        	}
        }
        else {
        	switch (i){
        		case (0) :
        			return 1 - x_r.x - x_r.y; 
        		case (1) :
        			return x_r.x;
        		case (2) : 
        			return x_r.y;
        	}	
        }
        return 0 ; // should not be reached
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        //std::cout << "[ShapeFunctions] evaluate gradient shape function " << i << '\n';
        vec2 g ;
        if (dim_ == 1){
        	switch (i){
        		case (0) : 
        			g.x = -1;
        			g.y = 0;
        			break;
        		case (1) : 
        			g.x = 1;
        			g.y = 0;
        			break;
        	}
        }
        else {
        	switch (i){
        		case (0) :
        			g.x = -1;
        			g.y = -1;
        			break;
        		case (1) :
        			g.x = 1;
        			g.y = 0;
        			break;
        		case (2) : 
        			g.x = 0;
        			g.y = 1;
        			break;
        	}	
        }
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),
        DenseMatrix& Ke )
    {
        //std::cout << "compute elementary matrix" << '\n';
        for(int i=0; i<3; ++i){
        
        	for(int j=0; j<3; ++j){
        	
        		double sum = 0;
        		
        		for(int pt =0; pt<quadrature.nb_points(); ++pt){
        		
        			vertex point = quadrature.point(pt);
        			double s1 = quadrature.weight(pt)*(*coefficient)(elt_mapping.transform(point));
        			DenseMatrix jacob = elt_mapping.jacobian_matrix(point);
        			vec2 gradphi_i = reference_functions.evaluate_grad(i, point);
        			vec2 gradphi_j = reference_functions.evaluate_grad(j, point);
        			vertex vec1 = ((jacob.invert_2x2()).transpose()).mult_2x2_2(gradphi_i);
        			vertex vec12 = ((jacob.invert_2x2()).transpose()).mult_2x2_2(gradphi_j);
        			double pdtscal = dot(vec1,vec12);
        			double detj = jacob.det_2x2();
        			sum = sum + (s1*pdtscal*detj);
        		}
        		
        		Ke.set(i, j, sum);
        	}
        }
        
        
        
        
        
    }

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
        //std::cout << "Ke -> K" << '\n';
                
        for(int i=0; i<Ke.height(); ++i){
        	for(int j=0; j<Ke.width(); ++j){
        		int i_global = M.get_triangle_vertex_index(t, i);
        		int j_global = M.get_triangle_vertex_index(t, j);
        		K.add(i_global, j_global, Ke.get(i, j));
        	}

        }
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        for(int i=0; i< reference_functions.nb_functions(); ++i){
        	
        		double pdt = 0;
        		
        		for(int pt =0; pt<quadrature.nb_points(); ++pt){
        		
        			vertex point = quadrature.point(pt);
        			double w = quadrature.weight(pt);
        			DenseMatrix jacob = elt_mapping.jacobian_matrix(point);
        			double phi_i = reference_functions.evaluate(i, point);
        			double detj = jacob.det_2x2();
        			pdt = pdt + w*phi_i*(*source)(elt_mapping.transform(point))*detj;
        		}
        		
        		Fe[i] = pdt;
        }
        
    }

    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (neumann condition)" << '\n';
        // TODO
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        if (border){
        	for(int ind=0; ind<Fe.size(); ++ind){
        		int i_global = M.get_edge_vertex_index(i, ind);
        		F[i_global] += Fe[ind];
       		}
        }
        else{
        	for(int ind=0; ind<Fe.size(); ++ind){
        		int i_global = M.get_triangle_vertex_index(i, ind);
        		F[i_global] += Fe[ind];
       		}


        }
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        //std::cout << "apply dirichlet boundary conditions" << '\n';
        
        std::vector<bool> processed_vertices(values.size(), false);
        double penalty_coefficient = 10000;
        for(int edge=0; edge<M.nb_edges(); edge++){
        	int edge_attribute = M.get_edge_attribute(edge);
        	if(attribute_is_dirichlet[edge_attribute]) {
        		for( int v=0; v<2; v++){
        			int vertex_index = M.get_edge_vertex_index(edge, v);
        			if( !processed_vertices[vertex_index]){
        				processed_vertices[vertex_index] = true;
        				K.add(vertex_index, vertex_index, penalty_coefficient);
        				F[vertex_index] += penalty_coefficient*values[vertex_index];
        			}
        		}
        	}
        }
        
    }

    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
    }

}
