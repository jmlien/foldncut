#include "datastructure_ss.h"

typedef Ss::Vertex_const_iterator Vertex_const_iterator ;
typedef Ss::Vertex::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;

template<class K>
void generate_perpendiculars(CGAL::Straight_skeleton_2<K> const& iss, CGAL::Straight_skeleton_2<K> const& ess, std::vector<stitchingGraph> const& stg, std::list<Halfedge_handle>& ppd)
{
	//*****init safety 
	//*****also perpendiculars for exterior skeleton
	//Loop through all the vertices 
	for(Vertex_const_iterator v = iss.vertices_begin(); v != iss.vertices_end(); ++v)
	{
		if(v->is_skeleton()) //if given vertex is a skeleton vertex,
		{
			//Loop through all the faces around vertex
			Halfedge_around_vertex_circulator e_root =  v->halfedge_around_vertex_begin();
			Halfedge_around_vertex_circulator e = v->halfedge_around_vertex_begin();
			std::cout <<"v: "<< v->point().x()<<" "<< v->point().y()<<std::endl;

			do{
				std::cout << "f: "<< (*e)->face()->id() << std::endl;
				Halfedge_handle edge_ss = (*e)->face()->halfedge();
				Halfedge_handle edge_ss_root = (*e)->face()->halfedge();
				do{
					//std::cout << edge_ss->vertex()->point().x() <<" "<< edge_ss->vertex()->point().y() << std::end;

					edge_ss = edge_ss->next();
				}while( edge_ss != edge_ss_root);
				//				Halfedge_handle first_ray;	//*****
				//				ppd.push_back(first_ray);
				//				perpendicular(iss, ess, stg, ppd);
				std::cout << std::endl;

				e++;
			}while( e != e_root);


		}
	}
}

template<class K>
void perpendicular(CGAL::Straight_skeleton_2<K> const& iss, CGAL::Straight_skeleton_2<K> const& ess, std::vector<stitchingGraph> const& stg, std::list<Halfedge_handle>& ppd)
{
	//***** Termination condition ----- add n?

	Halfedge_handle last_ray = ppd.pop_back();

}

//void mark_crease();