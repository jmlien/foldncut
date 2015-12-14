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
							std::cout << "skeleton vertex: "<< skeleton_v.x() <<" " << skeleton_v.y()<<std::endl;
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
	//Starts from the perpendicular hits the cutedge with previous skeleton face id information
	std::cout << "input fid: "<<prev_fid<<std::endl;
	//Argument Safety
	if(ppd.empty())	{std::cout << "Cannot compute perpendiculars recursively. Please check the first rays of give model." <<std::endl; return; }

	Segment last_seg = ppd.back();
	Point last_begin = last_seg.vertex(0);
	Point last_end = last_seg.vertex(1);	//This point should be on the cutedge.

	int offset = iss.size_of_faces();

	//Termination condition
	if(last_end.x() < minX-PAPER_THRESHOLD || last_end.x() > maxX+PAPER_THRESHOLD) return;	//If the new point is outside of the paper, terminates.
	if(last_end.y() < minY-PAPER_THRESHOLD || last_end.y() > maxY+PAPER_THRESHOLD) return;
	if(level > 100) return;																	//If the recursive level is greater than 50, terminates.
	if(prev_fid == -1) { std::cout << "wrong face id" << std::endl; return;}				//Wrong previous face id

	BridgingGraph connect_info;								

	int current_fid = -1;	
	bool hits_skeleton_edge = false;

	Line incident_ray = Line(last_begin, last_end);
	Line normal;
	Line reflected_ray;
	Point incidence_p(-1, -1);

	Halfedge_const_handle edge_around_face_root, edge_around_face;

	if( prev_fid >= offset)	
	{
		//Previous perpendicular hits from the exterior skeleton
		//New perpendicular will be add inside interior skeleton. 

		//Obtaining the bridging information between two faces
		search_bridgegraph<K>(bg, prev_fid, offset, connect_info);
		if(connect_info.fid_ess != prev_fid) {std::cout <<"Cannot find matching information in bridging graph." << std::endl; return;}
		current_fid = connect_info.fid_iss;		//New skeleton face id to advance

		//Using face id, find the skeleton face structure and loop throught all the edges
		find_skeleton_face(iss, current_fid, 0, edge_around_face_root);
		edge_around_face = Halfedge_const_handle(edge_around_face_root);

	}
	else		
	{
		//Previous perpendicular hits from the interior skeleton
		//New perpendicular will be add inside exterior skeleton. 

		//Obtaining the bridging information between two faces
		search_bridgegraph<K>(bg, prev_fid, 0, connect_info);
		if(connect_info.fid_iss != prev_fid) {std::cout <<"Cannot find matching information in bridging graph." << std::endl; return;}
		current_fid = connect_info.fid_ess;		//New skeleton face id to advance

		//Using face id, find the skeleton face structure and loop throught all the edges
		find_skeleton_face(ess, current_fid, offset, edge_around_face_root);
		edge_around_face = Halfedge_const_handle(edge_around_face_root);
	}
	std::cout << "middle face:" << current_fid <<std::endl;
	do{
		//Find the skeleton edge which intersect with incident ray
		if(edge_around_face->is_bisector())
		{
			Segment edge_ss(Point(edge_around_face->vertex()->point().x(), edge_around_face->vertex()->point().y()),
				Point(edge_around_face->opposite()->vertex()->point().x(), edge_around_face->opposite()->vertex()->point().y()));

			//Find the intersection of skeleton edge and incident ray
			CGAL::Object res_intersect = CGAL::intersection(edge_ss, incident_ray);
			if(assign(incidence_p, res_intersect))
			{
				//If there is intersection, keep the segment connect from cutedge to skeleton edge
				normal = Line(edge_ss.vertex(0), edge_ss.vertex(1));
				Segment perpendicular_first(last_end, incidence_p);
				ppd.push_back(perpendicular_first);
				hits_skeleton_edge = true;
				break;
			}
		}
		edge_around_face = edge_around_face->next();
	}while( edge_around_face != edge_around_face_root);

	//If there is no intersection, finish the perpendicular function
	if(!hits_skeleton_edge)	return; 

	//Compute the segment connect from skeleton edge to cutedge
	Point proj_last_end = normal.projection(last_end);
	Point reflect_last_end= last_end + 2*(proj_last_end-last_end);
	Segment perpendicular_second(incidence_p,reflect_last_end);
	ppd.push_back(Segment(incidence_p,reflect_last_end));

	//Obtain current face id
	edge_around_face = edge_around_face->opposite();
	if(edge_around_face->face() != NULL)
	{
		if( prev_fid >= offset)
			current_fid = edge_around_face->face()->id();
		else
			current_fid = edge_around_face->face()->id()+offset;
	}
	else current_fid = -1;
	if(current_fid == -1) return;
	std::cout << "middle2 face: "<<current_fid<<std::endl;

	//Make sure if this segment hits the cutedge. If not, make a segment hits straight skeleton
	bool hits_skeleton = false;
	BridgingGraph next_b;

	Point hit_p;
	Segment ppd_candidate = Segment(perpendicular_second);
	while(1)
	{
		//Obtain the cutedge
		if( current_fid >= offset)	search_bridgegraph<K>(bg, current_fid, offset, next_b);
		else						search_bridgegraph<K>(bg, current_fid, 0, next_b);

		//Test if given candidate hits the cut edge
		CGAL::Object res_hit = CGAL::intersection(next_b.cut_edge, ppd_candidate);
		if(assign(hit_p, res_hit))	break;	//if this candidate hits cut edge, move to next recursive step.
		else								//if this candidate doesn't hit cut edge, find other skeleton edge which intersect with it
		{
			if(ppd_candidate.vertex(1).x() < minX-PAPER_THRESHOLD || ppd_candidate.vertex(1).x() > maxX+PAPER_THRESHOLD) return;	//If the new point is outside of the paper, terminates.
			if(ppd_candidate.vertex(1).y() < minY-PAPER_THRESHOLD || ppd_candidate.vertex(1).y() > maxY+PAPER_THRESHOLD) return;
			
			//Find the new incident ray using last data
			Line new_incident_ray(ppd_candidate.vertex(0), ppd_candidate.vertex(1));
			
			//Find the face using last face data
			Halfedge_const_handle e_root, e;
			if( current_fid >= offset)	
				find_skeleton_face(ess, current_fid, offset, e_root);
			else
				find_skeleton_face(iss, current_fid, 0, e_root);

			e = Halfedge_const_handle(e_root);

			Point new_inc_p;		
			do{
				//Loop through the skeleton edges
				if(e->is_bisector())
				{
					Segment edge_ss(Point(e->vertex()->point().x(), e->vertex()->point().y()), Point(e->opposite()->vertex()->point().x(), e->opposite()->vertex()->point().y()));
					
					//Test if the skeleton edge has intersection with candidate
					CGAL::Object res_hit_new = CGAL::intersection(edge_ss, ppd_candidate);
					if(assign(new_inc_p, res_hit_new))	//if it has intersection
					{
						if(abs((new_inc_p-ppd_candidate.vertex(0)).squared_length())>1e-16)	//and if this point is not on the skeleton edge already crossed, 
						{	
							std::cout << "find new incident point" << new_inc_p.x() << " " << new_inc_p.y()<< std::endl;
							std::cout << "original point"<<ppd_candidate.vertex(0).x() << " " << ppd_candidate.vertex(0).y() << std::endl;
							ppd.pop_back();																		//remove candidate
							ppd_candidate = Segment(ppd_candidate.vertex(0), new_inc_p);
							ppd.push_back(ppd_candidate);														//add modified candidate
							normal = Line(edge_ss.vertex(0), edge_ss.vertex(1));
							std::cout << "line:" <<edge_ss.vertex(0).x()<<" "<<edge_ss.vertex(0).y()<<" "<<edge_ss.vertex(1).x()<<edge_ss.vertex(1).y()<<std::endl;
							hits_skeleton_edge = true;

							//Obtain current face id
							e = e->opposite();
							if(e->face() != NULL)
							{
								if( current_fid >= offset)
									current_fid = e->face()->id()+offset;
								else
									current_fid = e->face()->id();
							}
							else current_fid = -1;
							std::cout <<"new face id:"<< current_fid << std::endl;
							if(current_fid == -1) return;
							

							//Compute the segment connect from skeleton edge to cutedge
							Point proj_p = normal.projection(ppd_candidate.vertex(0));
							std::cout << "proj p:"<<proj_p.x() <<" " <<proj_p.y()<< std::endl;
							Point reflect_p= last_end + 2*(proj_p-ppd_candidate.vertex(0));
							std::cout << "reflect_p p:"<<reflect_p.x() <<" " <<reflect_p.y()<< std::endl;
							ppd_candidate = Segment(new_inc_p, reflect_p);
							ppd.push_back(ppd_candidate);

							break;
						}
					}
				}
				e = e->next();
			}while(e != e_root);

			break;
		}
	}

	//Computing perpendiculars recursively
	perpendicular(iss, ess, bg, ppd, current_fid, ++level);
}

template <class K>
void reflect(Segment const& incident_seg, Point& incidence_p, Halfedge_const_handle e, std::list<Segment> ppd, bool hits_skeleton_edge)
{
	Line incident_ray = Line(incident_seg.vertex(0), incident_seg.vertex(1));
	Halfedge_const_handle e_root = Halfedge_const_handle(e);
	Line normal;
	do{
		//Find the skeleton edge which intersect with incident ray
		if(e->is_bisector())
		{
			Segment edge_ss(Point(e->vertex()->point().x(), e->vertex()->point().y()),
				Point(e->opposite()->vertex()->point().x(), e->opposite()->vertex()->point().y()));

			//Find the intersection of skeleton edge and incident ray
			CGAL::Object res_intersect = CGAL::intersection(edge_ss, incident_ray);
			if(assign(incidence_p, res_intersect))
			{
				//If there is intersection, keep the segment connect from cutedge to skeleton edge
				normal = Line(edge_ss.vertex(0), edge_ss.vertex(1));
				ppd.push_back(Segment(incident_seg.vertex(1), incidence_p));
				hits_skeleton_edge = true;
				break;
			}
		}
		e = e->next();
	}while( e != e_root);
}

template<class K>
void perpendicular_new(CGAL::Straight_skeleton_2<K> const& iss, CGAL::Straight_skeleton_2<K> const& ess, std::vector<BridgingGraph> const& bg, std::list<Segment>& ppd, int prev_fid, bool hit_cutedge, int level)	
{
	//Argument Safety
	if(ppd.empty())	{std::cout << "Cannot compute perpendiculars recursively. Please check the first rays of give model." <<std::endl; return; }

	Segment last_seg = ppd.back();
	Point last_begin = last_seg.vertex(0);
	Point last_end = last_seg.vertex(1);	//This point should be on the cutedge.

	int offset = iss.size_of_faces();

	//Termination condition
	if(last_end.x() < minX-PAPER_THRESHOLD || last_end.x() > maxX+PAPER_THRESHOLD) return;	//If the new point is outside of the paper, terminates.
	if(last_end.y() < minY-PAPER_THRESHOLD || last_end.y() > maxY+PAPER_THRESHOLD) return;
	if(level > 100) return;																	//If the recursive level is greater than 50, terminates.
	if(prev_fid == -1) { std::cout << "wrong face id" << std::endl; return;}				//Wrong previous face id

	int current_fid = -1;	
	if(hit_cutedge)
	{
		BridgingGraph connect_info;								
		bool hits_skeleton_edge = false;

		Line incident_ray = Line(last_begin, last_end);
		Point incidence_p(-1, -1);

		Halfedge_const_handle edge_around_face_root, edge_around_face;

		if( prev_fid >= offset)	
		{
			//Previous perpendicular hits from the exterior skeleton
			//New perpendicular will be add inside interior skeleton. 

			//Obtaining the bridging information between two faces
			search_bridgegraph<K>(bg, prev_fid, offset, connect_info);
			if(connect_info.fid_ess != prev_fid) {std::cout <<"Cannot find matching information in bridging graph." << std::endl; return;}
			current_fid = connect_info.fid_iss;		//New skeleton face id to advance

			//Using face id, find the skeleton face structure and loop throught all the edges
			find_skeleton_face(iss, current_fid, 0, edge_around_face_root);
			edge_around_face = Halfedge_const_handle(edge_around_face_root);

		}
		else		
		{
			//Previous perpendicular hits from the interior skeleton
			//New perpendicular will be add inside exterior skeleton. 

			//Obtaining the bridging information between two faces
			search_bridgegraph<K>(bg, prev_fid, 0, connect_info);
			if(connect_info.fid_iss != prev_fid) {std::cout <<"Cannot find matching information in bridging graph." << std::endl; return;}
			current_fid = connect_info.fid_ess;		//New skeleton face id to advance

			//Using face id, find the skeleton face structure and loop throught all the edges
			find_skeleton_face(ess, current_fid, offset, edge_around_face_root);
			edge_around_face = Halfedge_const_handle(edge_around_face_root);
		}

		do{
			//Find the skeleton edge which intersect with incident ray
			if(edge_around_face->is_bisector())
			{
				Segment edge_ss(Point(edge_around_face->vertex()->point().x(), edge_around_face->vertex()->point().y()),
					Point(edge_around_face->opposite()->vertex()->point().x(), edge_around_face->opposite()->vertex()->point().y()));

				//Find the intersection of skeleton edge and incident ray
				CGAL::Object res_intersect = CGAL::intersection(edge_ss, incident_ray);
				if(assign(incidence_p, res_intersect))
				{
					//If there is intersection, keep the segment connect from cutedge to skeleton edge
					Segment perpendicular_first(last_end, incidence_p);
					ppd.push_back(perpendicular_first);
					hits_skeleton_edge = true;
					break;
				}
			}
			edge_around_face = edge_around_face->next();
		}while( edge_around_face != edge_around_face_root);

		if(!hits_skeleton_edge)	return; 

	}//HIT_CUT_EDGE
	else
	{
		Halfedge_const_handle e_root, e;

		if( prev_fid >= offset)	
			find_skeleton_face(ess, current_fid, offset, e_root);
		else
			find_skeleton_face(iss, current_fid, 0, e_root);

		e = Halfedge_const_handle(e_root);


		Point proj_last_end = normal.projection(last_end);
		Point reflect_last_end= last_end + 2*(proj_last_end-last_end);
		Segment perpendicular_second(incidence_p,reflect_last_end);
		ppd.push_back(Segment(incidence_p,reflect_last_end));

	}//HIT_SKELETON_EDGE
	//Computing perpendiculars recursively
	perpendicular_new(iss, ess, bg, ppd, current_fid, true, ++level);
}

#endif