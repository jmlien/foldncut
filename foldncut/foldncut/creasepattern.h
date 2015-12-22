#ifndef SS_CREASE_H
#define SS_CREASE_H
#include "datastructure_ss.h"

#include <CGAL/intersections.h>

//For loop
typedef Ss::Vertex_const_iterator Vertex_const_iterator;
typedef Ss::Vertex::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;
typedef std::list<Segment>::iterator SEG_I;
typedef std::list<Perpendiculars>::iterator PER_I;

double ratio;	//diagonal length

template<class K>
void generate_perpendiculars(CGAL::Straight_skeleton_2<K> const& iss, CGAL::Straight_skeleton_2<K> const& ess, std::vector<BridgingGraph> const& bg, std::list<Perpendiculars>& ppd)
{
	//Argument safty check
	std::cout << "Genetrating perpendiculars.." << std::endl;
	if(iss.vertices_begin() == iss.vertices_end()) {std::cout << "There is no interior skeleton. Please check interior skeleton again"<<std::endl; return;}
	if(ess.vertices_begin() == ess.vertices_end()) {std::cout << "There is no exterior skeleton. Please check exterior skeleton again"<<std::endl; return;}
	if(bg.empty()) { std::cout << "There is no bridging graph. Please check bridging graph again." << std::endl; return; }
	if(!ppd.empty()) ppd.clear();

	ratio = sqrt((maxX-minX)*(maxX-minX)+(maxY-minY)*(maxY-minY));

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
							Perpendiculars first_ppd;
							first_ppd.seg = Segment(skeleton_v, intersect_p);
							first_ppd.level = 0;

							ppd.push_back(first_ppd);
							//std::cout << "*******skeleton vertex: "<< skeleton_v.x() <<" " << skeleton_v.y()<<std::endl;
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

	//FOR EXTERIOR SKELETON
	int offset = iss.size_of_faces();
	for(Vertex_const_iterator v = ess.vertices_begin(); v != ess.vertices_end(); ++v)
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
							if(skeleton_v.x() < minX - PAPER_THRESHOLD || skeleton_v.x() > maxX + PAPER_THRESHOLD) break;
							if(skeleton_v.y() < minY - PAPER_THRESHOLD || skeleton_v.y() > maxY + PAPER_THRESHOLD) break;
							//std::cout << "*******skeleton vertex: "<< skeleton_v.x() <<" " << skeleton_v.y()<<std::endl;
							Perpendiculars first_ppd;
							first_ppd.seg =Segment(skeleton_v, intersect_p);
							first_ppd.level = 0;

							ppd.push_back(first_ppd);
							//std::cout << "face id: "<< edge_sf->face()->id() +offset << std::endl;
							perpendicular(iss, ess, bg, ppd,edge_sf->face()->id() + offset,0);
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
}

template<class K>
void perpendicular(CGAL::Straight_skeleton_2<K> const& iss, CGAL::Straight_skeleton_2<K> const& ess, std::vector<BridgingGraph> const& bg, std::list<Perpendiculars>& ppd, int prev_fid, int level)	
{
	//Starts from the perpendicular hits the cutedge with previous skeleton face id information
	//std::cout << "input fid: "<<prev_fid<<std::endl;
	//Argument Safety
	if(ppd.empty())	{std::cout << "Cannot compute perpendiculars recursively. Please check the first rays of give model." <<std::endl; return; }

	Segment last_seg = ppd.back().seg;
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
	//std::cout << "middle face:" << current_fid <<std::endl;
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

				Perpendiculars real_ppd;
				real_ppd.seg = Segment(last_end, incidence_p);
				real_ppd.level = level++;

				ppd.push_back(real_ppd);

				hits_skeleton_edge = true;
				break;
			}else
			{
				if(CGAL::do_intersect(edge_ss, incident_ray)) //*******************?!?!?!?!?!?
				{
					CGAL::Object res_intersect = CGAL::intersection(Line(edge_ss.vertex(0), edge_ss.vertex(1)), incident_ray);
					if(assign(incidence_p, res_intersect))
					{
						//If there is intersection, keep the segment connect from cutedge to skeleton edge
						normal = Line(edge_ss.vertex(0), edge_ss.vertex(1));
						Segment perpendicular_first(last_end, incidence_p);

						Perpendiculars real_ppd;
						real_ppd.seg = Segment(last_end, incidence_p);
						real_ppd.level = level++;

						ppd.push_back(real_ppd);
						hits_skeleton_edge = true;
						break;
					}
				}
			}
		}
		edge_around_face = edge_around_face->next();
	}while( edge_around_face != edge_around_face_root);

	//If there is no intersection, finish the perpendicular function
	if(!hits_skeleton_edge)	return; 

	//Compute the segment connect from skeleton edge to cutedge
	Point proj_last_end = normal.projection(last_end);
	Point reflect_last_end= last_end + 2*(proj_last_end-last_end);

	Perpendiculars unknown_ppd;
	unknown_ppd.seg = Segment(incidence_p,reflect_last_end);
	unknown_ppd.level = level++;

	ppd.push_back(unknown_ppd);

	//Obtain current face id
	edge_around_face = edge_around_face->opposite();
	if(edge_around_face->face() != NULL)
	{
		if( prev_fid >= offset)
			current_fid = edge_around_face->face()->id();
		else
			current_fid = edge_around_face->face()->id()+offset;
	}
	//std::cout << "middle2 face: "<<current_fid<<std::endl;

	//Make sure if this segment hits the cutedge. If not, make a segment hits straight skeleton
	bool hits_next_skeleton = false;
	BridgingGraph next_b;

	Point hit_p;
	Segment ppd_candidate = Segment(unknown_ppd.seg);
	while(1)
	{
		//Obtain the cutedge
		if( current_fid >= offset)	
		{
			bool found = search_bridgegraph<K>(bg, current_fid, offset, next_b);
			if(!found) return;
		}
		else						
		{
			search_bridgegraph<K>(bg, current_fid, 0, next_b);
		}
		//std::cout <<"while  current face id : "<< current_fid<<std::endl;
		//std::cout <<"current ppd:"<<ppd_candidate.vertex(0).x() << " " <<ppd_candidate.vertex(0).y() << " "<< ppd_candidate.vertex(1).x() << " " <<ppd_candidate.vertex(1).y()<<std::endl;

		//Find the new incident ray using last data
		Line new_incident_ray(ppd_candidate.vertex(0), ppd_candidate.vertex(1));

		//Test if given candidate hits the cut edge
		CGAL::Object res_hit = CGAL::intersection(next_b.cut_edge, new_incident_ray);
		if(assign(hit_p, res_hit))//if this candidate hits cut edge, move to next recursive step.
		{	
			if(abs((hit_p-ppd_candidate.vertex(1)).squared_length())>1e-16)
			{
				level--;
				ppd.pop_back();

				Perpendiculars real_ppd;
				real_ppd.seg = Segment(ppd_candidate.vertex(0), hit_p);
				real_ppd.level = level++;
				ppd.push_back(real_ppd);
			}
			break;
		}	
		else								//if this candidate doesn't hit cut edge, find other skeleton edge which intersect with it
		{
			if(ppd_candidate.vertex(0).x() < minX-PAPER_THRESHOLD || ppd_candidate.vertex(0).x() > maxX+PAPER_THRESHOLD) {ppd.pop_back(); return;}	//If the new point is outside of the paper, terminates.
			if(ppd_candidate.vertex(0).y() < minY-PAPER_THRESHOLD || ppd_candidate.vertex(0).y() > maxY+PAPER_THRESHOLD) {ppd.pop_back(); return;}

			if(ppd_candidate.vertex(1).x() < minX-PAPER_THRESHOLD || ppd_candidate.vertex(1).x() > maxX+PAPER_THRESHOLD) return;	//If the new point is outside of the paper, terminates.
			if(ppd_candidate.vertex(1).y() < minY-PAPER_THRESHOLD || ppd_candidate.vertex(1).y() > maxY+PAPER_THRESHOLD) return;

			//Find the face using last face data
			Halfedge_const_handle e_root, e;
			if( current_fid >= offset)	
				find_skeleton_face(ess, current_fid, offset, e_root);
			else
				find_skeleton_face(iss, current_fid, 0, e_root);

			e = Halfedge_const_handle(e_root);

			Point new_inc_p;	hits_next_skeleton=false;	
			do{
				//Loop through the skeleton edges
				if(e->is_bisector())
				{
					Segment edge_ss(Point(e->vertex()->point().x(), e->vertex()->point().y()), Point(e->opposite()->vertex()->point().x(), e->opposite()->vertex()->point().y()));

					//Test if the skeleton edge has intersection with candidate
					if(!CGAL::do_intersect(edge_ss, new_incident_ray)) { e = e->next(); continue;}	//**************?!?!?!?!?!?

					CGAL::Object res_hit_new = CGAL::intersection(edge_ss, new_incident_ray);
					if(assign(new_inc_p, res_hit_new))	//if it has intersection
					{
						if(abs((new_inc_p-ppd_candidate.vertex(0)).squared_length())>1e-16)	//and if this point is not on the skeleton edge already crossed, 
						{	
							//std::cout << "skeleton edge: " << edge_ss.vertex(0).x() <<" " <<edge_ss.vertex(0).y() <<" "<<edge_ss.vertex(1).x()<<" "<< edge_ss.vertex(1).y()<<std::endl;
							//std::cout << "line: "<<ppd_candidate.vertex(0).x()<<" "<<ppd_candidate.vertex(0).y()<<" "<<ppd_candidate.vertex(1).x()<<" "<<ppd_candidate.vertex(1).y()<<std::endl;
							//std::cout << "find new incident point" << new_inc_p.x() << " " << new_inc_p.y()<< std::endl;
							//std::cout << "original point"<<ppd_candidate.vertex(0).x() << " " << ppd_candidate.vertex(0).y() << std::endl;
							//std::cout << ppd_candidate.vertex(1).x() << " " << ppd_candidate.vertex(1).y() << std::endl;
							level--;
							ppd.pop_back();																		//remove candidate
							ppd_candidate = Segment(ppd_candidate.vertex(0), new_inc_p);

							Perpendiculars imaginary_ppd;
							imaginary_ppd.seg = ppd_candidate;
							imaginary_ppd.level = level++;

							ppd.push_back(imaginary_ppd);														//add modified candidate
							normal = Line(edge_ss.vertex(0), edge_ss.vertex(1));
							//std::cout << "line:" <<edge_ss.vertex(0).x()<<" "<<edge_ss.vertex(0).y()<<" "<<edge_ss.vertex(1).x()<<edge_ss.vertex(1).y()<<std::endl;
							hits_next_skeleton = true;

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
							//std::cout <<"new face id:"<< current_fid << std::endl;
							if(current_fid == -1) return;

							//Compute the segment connect from skeleton edge to cutedge
							Point proj_p = normal.projection(ppd_candidate.vertex(0));
							//std::cout << "proj p:"<<proj_p.x() <<" " <<proj_p.y()<< std::endl;
							Point reflect_p= ppd_candidate.vertex(0) + 2*(proj_p-ppd_candidate.vertex(0));
							//std::cout << "reflect_p p:"<<reflect_p.x() <<" " <<reflect_p.y()<< std::endl;
							ppd_candidate = Segment(new_inc_p, reflect_p);

							Perpendiculars imaginary_ppd_2;
							imaginary_ppd_2.seg = ppd_candidate;
							imaginary_ppd_2.level = level++;

							ppd.push_back(imaginary_ppd_2);

							break;
						}//else std::cout<< "already crossed"<<std::endl;
					}
					else //TEST IF SKELETON EDGE INTERSECTION WITH CANDIDATES
					{
						if(CGAL::do_intersect(edge_ss, new_incident_ray)) //*************?!?!!?!?!?!!?!!?!?
						{
							//std::cout<<"!";
							CGAL::Object res_hit_new = CGAL::intersection(Line(edge_ss.vertex(0), edge_ss.vertex(1)), new_incident_ray);
							if(assign(new_inc_p, res_hit_new)){
								if(abs((new_inc_p-ppd_candidate.vertex(0)).squared_length())>1e-16)	//and if this point is not on the skeleton edge already crossed, 
								{	
									level--;
									ppd.pop_back();																		//remove candidate
									ppd_candidate = Segment(ppd_candidate.vertex(0), new_inc_p);

									Perpendiculars imaginary_ppd;
									imaginary_ppd.seg = ppd_candidate;
									imaginary_ppd.level = level++;

									ppd.push_back(imaginary_ppd);														//add modified candidate
									normal = Line(edge_ss.vertex(0), edge_ss.vertex(1));
									hits_next_skeleton = true;

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
									if(current_fid == -1) return;

									//Compute the segment connect from skeleton edge to cutedge
									Point proj_p = normal.projection(ppd_candidate.vertex(0));
									Point reflect_p= ppd_candidate.vertex(0) + 2*(proj_p-ppd_candidate.vertex(0));
									ppd_candidate = Segment(new_inc_p, reflect_p);

									Perpendiculars imaginary_ppd_2;
									imaginary_ppd_2.seg = ppd_candidate;
									imaginary_ppd_2.level = level++;

									ppd.push_back(imaginary_ppd_2);

									break;
								}
							}
						}
					}
				}
				e = e->next();
			}while(e != e_root);	//LOOP_THROUGH_CURRENT_FACE
			if(!hits_next_skeleton) return;
		}	//CANDIDATE DOESN'T HIT CUT EDGE
	}
	//Computing perpendiculars recursively
	perpendicular(iss, ess, bg, ppd, current_fid, level);
}

template <class K>
void deduplicate_perpendiculars(std::list<Perpendiculars>& ppd)
{

	//Erase if perpendiculars are too closed. This is due to numerical error from straight skeleton
	std::cout << "dedpulicate perpendiculars..."<<std::endl;

	for(PER_I i = ppd.begin() ; i != ppd.end(); ++i)
	{
		PER_I j = i; ++j;		
		for( ; j != ppd.end(); ++j)
		{
			if(idential_segment_test(i->seg, j->seg, ratio)) 
			{
				if(j->level < i->level) 
				{
					i = ppd.erase(i);
					if(i==ppd.end()) break;
					else
					{
						j = i; ++j;
						continue;
					}
				}
				j = ppd.erase(j);
				if(j== ppd.end()) break;
				else continue;
			}
		}
	}


	//Erase the broken perpendiculars pieces
	PER_I j;
	for(PER_I i = ppd.begin() ; i != ppd.end();++i)
	{
		j = i; ++j;
		if(j == ppd.end()) continue;
		while(j->level - i->level > 1)
		{
			j = ppd.erase(j);
			if(j == ppd.end()) break;
		}
	}
}

template<class K>
void MountainValley(CGAL::Straight_skeleton_2<K> const& iss,CGAL::Straight_skeleton_2<K> const& ess, std::vector<BridgingGraph> const &bg, std::list<Perpendiculars> const& ppd, std::list<Segment>& mt, std::list<Segment>& vl)
{
	//**** can remove polygon
	std::cout << "Deciding mountain and valley..."<<std::endl;

	Segment in_edge;
	Segment out_edge;

	//LOOP INERIOR SKELETON
	for(Vertex_const_iterator v = iss.vertices_begin(); v != iss.vertices_end(); ++v)
	{
		//FOR VERTEX ON BOUNDARY
		if(v->is_contour()) //if given vertex is a contour vertex,
		{
			Vertex_const_iterator prev_v = v->prev_link;
			Vertex_const_iterator next_v = v->next_link;

			if(next_v->is_skeleton())
			{
				next_v =  iss.vertices_begin();
			}

			if(v == iss.vertices_begin())
			{
				prev_v = iss.vertices_end();
				do{
					--prev_v;
				}while(prev_v->is_skeleton());
			}

			in_edge = Segment(prev_v->point(), v->point());
			out_edge = Segment(v->point(),next_v->point());

			//MOUNTAIN
			if(CGAL::left_turn(in_edge.vertex(0), in_edge.vertex(1), out_edge.vertex(1)))	//If this vertex is convex vertex
			{
				//Find the skeleton edge around this vertex
				Halfedge_around_vertex_const_circulator e_root = v->halfedge_around_vertex_begin();
				Halfedge_around_vertex_const_circulator e = v->halfedge_around_vertex_begin();
				do	
				{
					if((*e)->is_bisector())
					{
						//Add skeleton edge as mountain
						mt.push_back(Segment((*e)->vertex()->point(), (*e)->opposite()->vertex()->point()));
					}
					++e;
				}while( e != e_root);
			}
			//VALLEY
			else		//this vertex is reflect vertex																	
			{
				//Find skeleton edge around this vertex
				Halfedge_around_vertex_const_circulator e_root = v->halfedge_around_vertex_begin();
				Halfedge_around_vertex_const_circulator e = v->halfedge_around_vertex_begin();
				do	
				{
					if((*e)->is_bisector())
					{
						//Add skeleton edge as valley
						vl.push_back(Segment((*e)->vertex()->point(), (*e)->opposite()->vertex()->point()));
					}
					++e;
				}while( e != e_root);
			}
		}
		//FOR STRAIGHT SKELETON VERTEX
		else
		{
			//Find the skeleton edge around this vertex
			Halfedge_around_vertex_const_circulator e_root = v->halfedge_around_vertex_begin();
			Halfedge_around_vertex_const_circulator e = v->halfedge_around_vertex_begin();
			do	
			{
				//If this sktraight skeleton edge consists of two straight skeleton vertex, add it as mountain.
				if((*e)->opposite()->vertex()->is_skeleton())
					mt.push_back(Segment((*e)->vertex()->point(), (*e)->opposite()->vertex()->point()));
				++e;
			}while( e != e_root);
		}
	}//LOOP INERIOR SKELETON

	//LOOP EXTERIOR SKELETON
	for(Vertex_const_iterator v = ess.vertices_begin(); v != ess.vertices_end(); ++v)
	{
		if(v->point().x() < minX-PAPER_THRESHOLD || v->point().x() > maxX+PAPER_THRESHOLD) continue;
		if(v->point().y() < minY-PAPER_THRESHOLD || v->point().y() > maxY+PAPER_THRESHOLD) continue;
	
		//FOR VERTEX ON BOUNDARY
		if(v->is_contour()) //if given vertex is a contour vertex,
		{
			Vertex_const_iterator prev_v = v->prev_link;
			Vertex_const_iterator next_v = v->next_link;
			
			in_edge = Segment(prev_v->point(), v->point());
			out_edge = Segment(v->point(), next_v->point());

			//MOUNTAIN
			if(CGAL::right_turn(in_edge.vertex(0), in_edge.vertex(1), out_edge.vertex(1)))	//If this vertex is convex vertex 
			{
				//Find the skeleton edge around this vertex
				Halfedge_around_vertex_const_circulator e_root = v->halfedge_around_vertex_begin();
				Halfedge_around_vertex_const_circulator e = v->halfedge_around_vertex_begin();
				do	
				{
					if((*e)->is_bisector())
					{
						//Add skeleton edge as mountain
						mt.push_back(Segment((*e)->vertex()->point(), (*e)->opposite()->vertex()->point()));
					}
					++e;
				}while( e != e_root);
			}
			//VALLEY
			else		//this vertex is reflect vertex																	
			{
				//Find skeleton edge around this vertex
				Halfedge_around_vertex_const_circulator e_root = v->halfedge_around_vertex_begin();
				Halfedge_around_vertex_const_circulator e = v->halfedge_around_vertex_begin();
				do	
				{
					if((*e)->is_bisector())
					{
						//Add skeleton edge as valley
						vl.push_back(Segment((*e)->vertex()->point(), (*e)->opposite()->vertex()->point()));
					}
					++e;
				}while( e != e_root);
			}
		}
		//FOR STRAIGHT SKELETON VERTEX
		//else
		//{
		//	//Find the skeleton edge around this vertex
		//	Halfedge_around_vertex_const_circulator e_root = v->halfedge_around_vertex_begin();
		//	Halfedge_around_vertex_const_circulator e = v->halfedge_around_vertex_begin();
		//	do	
		//	{
		//		//If this sktraight skeleton edge consists of two straight skeleton vertex, add it as mountain.
		//		if((*e)->opposite()->vertex()->is_skeleton())
		//			mt.push_back(Segment((*e)->vertex()->point(), (*e)->opposite()->vertex()->point()));
		//		++e;
		//	}while( e != e_root);
		//}
	}//LOOP EXTERIOR SKELETON
}
#endif