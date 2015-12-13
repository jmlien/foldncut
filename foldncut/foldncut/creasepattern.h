#ifndef SS_CREASE_H
#define SS_CREASE_H
#include "datastructure_ss.h"

#include <CGAL/intersections.h>

//For loop
typedef Ss::Vertex_const_iterator Vertex_const_iterator;
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

						Line cut_edge = Line(cut_edge_v1, cut_edge_v2);
						Line first_ray = cut_edge.perpendicular(skeleton_v);

						//Find the intersection point between cut edge and first ray starts from skeleton vertex. 
						CGAL::Object res_intersect = CGAL::intersection(Segment(cut_edge_v1, cut_edge_v2), first_ray);							//***** can be modifed to segment & line intersection
						Point intersect_p;
						if(assign(intersect_p, res_intersect))
						{
							ppd.push_back(Segment(skeleton_v, intersect_p));
							perpendicular(iss, ess, bg, ppd,edge_sf->face()->id(),0);
							break;
						}
						else ;//There is something wrong in computation. That causes intersection goes wrong
					}
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

	Segment last_seg = ppd.back();
	Point last_begin = last_seg.vertex(0);
	Point last_end = last_seg.vertex(1);	//This point should be on the cutedge.

	//Termination condition
	if(last_end.x() < minX-PAPER_THRESHOLD || last_end.x() > maxX+PAPER_THRESHOLD) return;	//If the new point is outside of the paper, terminates.
	if(last_end.y() < minY-PAPER_THRESHOLD || last_end.y() > maxY+PAPER_THRESHOLD) return;
	if(level > 50) return;																	//If the recursive level is greater than 50, terminates.
	if(prev_fid == -1) return;																//Wrong previous face id

	//Reflect ray
	BridgingGraph connect_info;								//Obtaining the bridging information between two faces
	int offset = iss.size_of_faces();
	if( prev_fid > offset)
		search_bridgegraph<K>(bg, prev_fid, offset, connect_info);
	else 
		search_bridgegraph<K>(bg, prev_fid, 0, connect_info);

	int current_fid = -1;
	bool incident_point_exists = false;

	if( prev_fid > offset)		//Reflects from exterior skeleton to interior skeleton 
	{
		Line incident_ray = Line(last_begin, last_end);
		Line normal;
		Line reflected_ray;
		Point incidence_p(-1, -1);

		if(connect_info.fid_ess != prev_fid) {std::cout <<"Cannot find matching information in bridging graph." << std::endl; return;}

		//Using face id, find the exterior skeleton face to advance
		current_fid = connect_info.fid_iss;
		Halfedge_const_handle edge_around_face_root, edge_around_face;
		find_skeleton_face(iss, current_fid, 0, edge_around_face_root);

		edge_around_face = Halfedge_const_handle(edge_around_face_root);
		do{
			if(edge_around_face->is_bisector())
			{
				Segment edge_ss(Point(edge_around_face->vertex()->point().x(), edge_around_face->vertex()->point().y()),
					Point(edge_around_face->opposite()->vertex()->point().x(), edge_around_face->opposite()->vertex()->point().y()));

				//Test intersection

				CGAL::Object res_intersect = CGAL::intersection(edge_ss, incident_ray);
				if(assign(incidence_p, res_intersect))
				{
					normal = Line(edge_ss.vertex(0), edge_ss.vertex(1));
					ppd.push_back(Segment(last_end, incidence_p));
					incident_point_exists = true;
					break;
				}
			}

			edge_around_face = edge_around_face->next();
		}while( edge_around_face != edge_around_face_root);

		if(!incident_point_exists) 
		{
			ppd.pop_back(); //***** need to implement 
			return;
		}
		
		edge_around_face = edge_around_face->opposite();
		current_fid = edge_around_face->face()->id();

		Point proj_last_end = normal.projection(last_end);
		Point reflect_last_end= last_end + 2*(proj_last_end-last_end);	
		ppd.push_back(Segment(incidence_p,reflect_last_end));
	}
	else		//Reflects from interior skeleton to exterior skeleton
	{
		Line incident_ray = Line(last_begin, last_end);
		Line normal;
		Line reflected_ray;
		Point incidence_p(-1, -1);
		
		if(connect_info.fid_iss != prev_fid) {std::cout <<"Cannot find matching information in bridging graph." << std::endl; return;}

		//Using face id, find the exterior skeleton face to advance
		Halfedge_const_handle edge_around_face_root, edge_around_face;
		current_fid = connect_info.fid_ess;
		find_skeleton_face(ess, current_fid, offset, edge_around_face_root);
	
		//Find the incident point and add incident segment to perpendiculars
		edge_around_face = Halfedge_const_handle(edge_around_face_root);
		do{
			if(edge_around_face->is_bisector())
			{
				Segment edge_ss(Point(edge_around_face->vertex()->point().x(), edge_around_face->vertex()->point().y()),
					Point(edge_around_face->opposite()->vertex()->point().x(), edge_around_face->opposite()->vertex()->point().y()));

				//Test intersection
				CGAL::Object res_intersect = CGAL::intersection(edge_ss, incident_ray);
				if(assign(incidence_p, res_intersect))
				{
					normal = Line(edge_ss.vertex(0), edge_ss.vertex(1));
					ppd.push_back(Segment(last_end, incidence_p));
					incident_point_exists = true;
					break;
				}
			}
			edge_around_face = edge_around_face->next();
		}while( edge_around_face != edge_around_face_root);

		if(!incident_point_exists)
		{
				ppd.pop_back();
				return;
		}
		
		edge_around_face = edge_around_face->opposite();
		if(edge_around_face->face() != NULL) current_fid = edge_around_face->face()->id()+offset;
		else current_fid = -1;

		//Compute reflected ray
		Point proj_last_end = normal.projection(last_end);
		Point reflect_last_end= last_end + 2*(proj_last_end-last_end);
		ppd.push_back(Segment(incidence_p,reflect_last_end));
	}

	//Computing perpendiculars recursively
	perpendicular(iss, ess, bg, ppd, current_fid, ++level);
}

#endif