#ifndef SS_DATASTRUCTURE_H
#define SS_DATASTRUCTURE_H

#include "globals.h"

#include <CGAL/circulator.h>

//For Bridging Graph
typedef Polygon_2::Edge_const_circulator Edge_const_circulator;
typedef Ss::Halfedge_const_iterator Halfedge_const_iterator ;
typedef Ss::Halfedge_handle Halfedge_handle;
typedef Ss::Halfedge_const_handle Halfedge_const_handle;
typedef Ss::Face_const_handle Face_const_handle;
typedef Ss::Face_handle Face_handle;

typedef struct BridgingGraph
{
	Segment cut_edge;		
	int fid_iss;
	int fid_ess;

}BRIDGINGGRAPH;

typedef std::vector<BridgingGraph>::const_iterator I;
typedef CGAL::Circulator_from_iterator<I>	Circulator;

typedef Ss::Face_const_iterator Face_const_iterator ;

bool is_almost_same_segment(Segment const& first, Segment const& second, double epsilon)
{
	int same_v = 0;
	Point fa(first.vertex(0)), fb(first.vertex(1));
	Point sa(second.vertex(0)), sb(second.vertex(1));

	double dist1 = abs((fa-sa).squared_length());
	double dist2 = abs((fa-sb).squared_length());
	double dist3 = abs((fb-sa).squared_length());
	double dist4 = abs((fb-sb).squared_length());

	if(dist1 < epsilon ) same_v++;
	if(dist2 < epsilon ) same_v++;
	if(dist3 < epsilon ) same_v++;
	if(dist4 < epsilon ) same_v++;

	if(same_v == 2)	return true;
	else return false;
}
bool idential_segment_test(Segment const& first, Segment const& second, double epsilon)
{
	int same_v = 0;
	Point fa(first.vertex(0)), fb(first.vertex(1));
	Point sa(second.vertex(0)), sb(second.vertex(1));

	double dist1 = abs((fa-sa).squared_length());
	double dist2 = abs((fa-sb).squared_length());
	double dist3 = abs((fb-sa).squared_length());
	double dist4 = abs((fb-sb).squared_length());

	if(dist1 < epsilon ) same_v++;
	if(dist2 < epsilon ) same_v++;
	if(dist3 < epsilon ) same_v++;
	if(dist4 < epsilon ) same_v++;

	if(same_v == 2)	
	{
		if( (dist1+dist4) < epsilon || (dist2 + dist3) < epsilon ) return true;
		else return false;
	}
	else return false;
}

template<class K>
void construct_bridging_graph(Polygon_2& poly, CGAL::Straight_skeleton_2<K> const& interior_ss, CGAL::Straight_skeleton_2<K> const& exterior_ss, std::vector<BridgingGraph>& bg)
{
	std::cout << "Constructing data structure.." << std::endl;

	if(poly.is_empty()) {std::cout << "Polygon is empty. Please check your model again."<< std::endl; exit(-1);}

	//Make circulator in the same order with polygon
	//This circulator contains information about cut edges and connected face to cut edges
	Edge_const_circulator root =  poly.edges_circulator();
	Edge_const_circulator e = poly.edges_circulator();
	do{
		BridgingGraph bgNode = BridgingGraph();
		bgNode.cut_edge = Segment(Point(e->vertex(0).x(), e->vertex(0).y()), Point(e->vertex(1).x(), e->vertex(1).y()));
		bgNode.fid_iss = find_skeleton_face_id(interior_ss, bgNode.cut_edge, 0);
		bgNode.fid_ess = find_skeleton_face_id(exterior_ss, bgNode.cut_edge, interior_ss.size_of_faces());

		//std::cout << "cut 1 :" << e->vertex(0).x() << std::endl;
		//std::cout << "bg iss :" << bgNode.fid_iss << " bg ess: " << bgNode.fid_ess << std::endl;

		if(bgNode.fid_iss == -1) {std::cout << "Cannot find match cutedges in straight skeleton structure" << std::endl; exit(-1);}
		if(bgNode.fid_ess == -1) {std::cout << "Cannot find match cutedges in straight skeleton structure" << std::endl; exit(-1);}
		bg.push_back(bgNode);
		++e;
	}while( e != root);

}

//This function finds the face id of exterior skeleton and interior skeleton contains the cut edge
template<class K>
int find_skeleton_face_id(CGAL::Straight_skeleton_2<K> const& ss, Segment& seg, int offset)
{
	//***** This function can be improved, becuase both polygon and straight skeleton structures are stored in particular order. 
	double epsilon = 1e-16;

	for(Halfedge_const_iterator h = ss.halfedges_begin(); h != ss.halfedges_end(); ++h)
	{	

		//Find if the edge is identical to halfedge in straight skeleton structure
		Segment seg_ss(Point(h->vertex()->point().x(), h->vertex()->point().y()),Point(h->opposite()->vertex()->point().x(), h->opposite()->vertex()->point().y()));

		if(idential_segment_test(seg, seg_ss, epsilon))
		{
			//Keep the id of skeleton face
			if(h->face() != NULL) return h->face()->id()+offset;
			else if(h->opposite()->face() != NULL) return h->opposite()->face()->id()+offset;
			else return -1;
		}
	}
	return -1;
}

//This function can help move on from cut edges to next skeleton faces(where reflected ray exists). 
template<class K>
void find_skeleton_face(CGAL::Straight_skeleton_2<K> const& ss, int const& fid, int const& offset, Halfedge_const_handle& edge_ss)
{
	Face_const_iterator f = ss.faces_begin();
	int n = fid - offset;
	while( n != 0 )
	{
		++f;
		--n;
	}
	edge_ss = f->halfedge();
}

template<class K>
void search_bridgegraph(std::vector<BridgingGraph> const& bg_list, Segment const& e, BridgingGraph& bg)
{
	//***** This function can be improved, because of special order in CGAL datastructure
	bool found = false;
	double epsilon = 1e-16;

	for(I b=bg_list.begin(); b != bg_list.end(); ++b)
	{
		if(idential_segment_test(e, b->cut_edge, epsilon)) 
		{ 
			bg.cut_edge=b->cut_edge;
			bg.fid_ess = b->fid_ess;
			bg.fid_iss = b->fid_iss;

			found = true;
			return;
		}
	}
	if(!found) std::cout << "Cannot find matched bridging graph!" << std::endl;
}

template<class K>
bool search_bridgegraph(std::vector<BridgingGraph> const& bg_list, int const& fid, int const& offset, BridgingGraph& bg)
{
	bool found = false;
	//SEARCH ON EXTERIOR SKELTON
	if( offset != 0)
	{
		for(I b=bg_list.begin(); b != bg_list.end(); ++b)
		{
			if(b->fid_ess == fid) 
			{ 
				bg.cut_edge = b->cut_edge;
				bg.fid_ess = b->fid_ess;
				bg.fid_iss = b->fid_iss;

				found = true;
				break;
			}
		}
	}
	//SEARCH ON INTERIOR SKELTON
	else
	{
		for(I b=bg_list.begin(); b != bg_list.end(); ++b)
		{
			if(b->fid_iss == fid) 
			{ 
				bg.cut_edge = b->cut_edge;
				bg.fid_ess = b->fid_ess;
				bg.fid_iss = b->fid_iss;

				found = true;
				break;
			}
		}
	}
	return found;
}
#endif