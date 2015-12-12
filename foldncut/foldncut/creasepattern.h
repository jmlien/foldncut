#ifndef SS_CREASE_H
#define SS_CREASE_H
#include "datastructure_ss.h"

#include <CGAL/intersections.h>

//For loop
typedef Ss::Vertex_const_iterator Vertex_const_iterator;
typedef Ss::Halfedge_const_handle Halfedge_const_handle;
typedef Ss::Vertex::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;

template<class K>
void generate_perpendiculars(CGAL::Straight_skeleton_2<K> const& iss, CGAL::Straight_skeleton_2<K> const& ess, std::vector<BridgingGraph> const& bg, std::list<Segment>& ppd)
{
	//Argument safty check
	std::cout << "Genetrating perpendiculars.." << std::endl;
	if(iss.vertices_begin() == iss.vertices_end()) {std::cout << "There is no interior skeleton. Please check interior skeleton again"<<std::endl; return;}
	if(ess.vertices_begin() == ess.vertices_end()) {std::cout << "There is no exterior skeleton. Please check exterior skeleton again"<<std::endl; return;}
	if(bg.empty()) { std::cout << "There is no bridging graph. Please check bridging graph again." << std::endl; return; }
	if(!ppd.empty()) ppd.clear();

	//FOR INTERIOR SKELETON
	//Loop through all the vertices 
	for(Vertex_const_iterator v = iss.vertices_begin(); v != iss.vertices_end(); ++v)
	{
		if(v->is_skeleton()) //if given vertex is a skeleton vertex,
		{
			std::cout << "degree: "<< v->degree() << std::endl;
			//Loop throguth all the skeleton faces around given vertex
			//Compute the first ray of each face and start computing perpendicualr recursively
			Halfedge_around_vertex_const_circulator e_root = v->halfedge_around_vertex_begin();
			Halfedge_around_vertex_const_circulator e = v->halfedge_around_vertex_begin();
			do
			{
				//Find the first edge which the first ray hits
				Halfedge_const_handle edge_sf = (*e)->face()->halfedge();
				Halfedge_const_handle edge_sf_root = (*e)->face()->halfedge();
				do{
					if(!(edge_sf->is_bisector()))	//if this edge is on the polygon
					{

						Point skeleton_v(v->point().x(), v->point().y());
						Point cut_edge_v1(edge_sf->vertex()->point().x(), edge_sf->vertex()->point().y());
						Point cut_edge_v2(edge_sf->opposite()->vertex()->point().x(), edge_sf->opposite()->vertex()->point().y());

						std::cout << "skeleton vertex: "<< skeleton_v.x()<< " "<<skeleton_v.y()<< std::endl;
						std::cout << "polygon edge: "<< edge_sf->vertex()->point().x() <<" "<< edge_sf->vertex()->point().y() <<" "<< edge_sf->opposite()->vertex()->point().x() << " " << edge_sf->opposite()->vertex()->point().y() << std::endl;

						Line cut_edge = Line(cut_edge_v1, cut_edge_v2);
						Line first_ray = cut_edge.perpendicular(skeleton_v);

						std::cout << "cut edge: "<< cut_edge.a()<<" "<<cut_edge.b()<<" "<<cut_edge.c()<<std::endl;
						std::cout << "first ray: "<< first_ray.a() <<" "<<first_ray.b()<<" "<<first_ray.c() << std::endl;

						//Find the intersection point between cut edge and first ray starts from skeleton vertex. 
						CGAL::Object res_intersect = CGAL::intersection(cut_edge, first_ray);							//***** can be modifed to segment & line intersection
						Point intersect_p;
						if(assign(intersect_p, res_intersect))
						{
							std::cout << "intersection: "<< intersect_p.x()<<" "<< intersect_p.y()<<std::endl;
							if(Segment(cut_edge_v1, cut_edge_v2).collinear_has_on(intersect_p))//If the first ray doesn't hit the cut edge, doe
								ppd.push_back(Segment(skeleton_v, intersect_p));
							else std::cout<< "intersection point isn't on the cut edge" << std::endl;
							break;
						}
						else ;//There is something wrong in computation. That causes intersection goes wrong
						std::cout << std::endl;
					}
					std::cout << std::endl;
					edge_sf = edge_sf->next();
				}while( edge_sf != edge_sf_root);
				++e;
			}while( e != e_root);
		}
	}

	//*****also perpendiculars for exterior skeleton
}

template<class K>
void perpendicular(CGAL::Straight_skeleton_2<K> const& iss, CGAL::Straight_skeleton_2<K> const& ess, std::vector<BridgingGraph> const& bg, std::list<Segment>& ppd, int prev_fid, int level)	
{
	//Argument Safety
	if(ppd.empty())	{std::cout << "Cannot compute perpendiculars recursively. Please check the first rays of give model." <<std::endl; return; }

	Segment last_ray = ppd.back();
	Point last_begin = last_ray.vertex(0);
	Point last_end = last_ray.vertex(1);	//This point should be on the cutedge.

	//Termination condition
	if(last_end.x() < minX || last_end.x() > maxX) return;	//If the new point is outside of the paper, terminates.
	if(last_end.y() < minY || last_end.y() > maxY) return;
	if(level > 50) return;										//If the recursive level is greater than 50, terminates.

	//Reflect ray
	BridgingGraph connect_info;								//Obtaining the bridging information between two faces
	search_bridgegraph(bg, prev_fid, connect_info);
	int offset = interior_ss.size_of_faces();

	Line incident_ray = Line(last_begin, last_end);
	Line reflected_ray;
	Angle incidence_angle;
	Point incidence_p;

	int current_fid = -1;

	if( prev_fid > offset)		//Reflects from exterior skeleton to interior skeleton 
	{
		if(connect_info.fid_ess != prev_fid) {std::cout <<"Cannot find matching information in bridging graph." << std::endl; return;}

		//Using face id, find the skeleton face to advance
		Halfedge_handle edge_around_face_root, edge_around_face;
		find_skeleton_face(ess, prev_fid, offset, edge_around_face_root);	//*****To edge const

		//Add incident segment
		edge_around_face = Halfedge_handle(edge_around_face_root);
		do{
			if(edge_around_face->is_bisector())
			{
				Segment edge_ss(Point(edge_around_face->vertex()->point().x(), edge_around->vertex()->point().y()),
					Point(edge_around_face->opposite()->vertex()->point().x(), edge_around->opposite()->vertex()->point().y()));

				//Test intersection
				CGAL::Object res_intersect = CGAL::intersection(edge_ss, incident_ray);
				Point incidence_p;
				if(assign(incidence_p, res_intersect))
				{
					//***** degenerate case?
					ppd.push_back(last_end, incidence_p);
					break;
				}
				else ;	//***** 
			}
		}while( edge_around_face != edge_around_face_root);

		//compute the angle of incidence and reflect the ray
		//incidence_angle = CGAL::angle(last_end, incidence_p, edge_ss

		//Add reflected segment

		//incidence_angle = CGAL::angle(last_begin, last_end, 
	}
	else										//Reflects from interior skeleton to exterior skeleton
	{
		if(connect_info.fid_iss != prev_fid) {std::cout <<"Cannot find matching information in bridging graph." << std::endl; return;}
	}

	//Computing perpendiculars recursively
	perpendicular(iss, ess, bg, ppd, current_fid, ++level);
}

#endif