#ifndef SS_DATASTRUCTURE_H
#define SS_DATASTRUCTURE_H

#include "globals.h"

#include <CGAL/circulator.h>

//For Stitched Graph

typedef Polygon_2::Edge_const_circulator Edge_const_circulator;
typedef Ss::Halfedge_const_iterator Halfedge_const_iterator ;
typedef Ss::Halfedge_handle Halfedge_handle;
typedef Ss::Halfedge_const_handle Halfedge_const_handle;

typedef struct BridgingGraph
{
	Segment cut_edge;		
	int fid_iss;
	int fid_ess;

}BRIDGINGGRAPH;
typedef std::vector<BridgingGraph>::const_iterator I;
typedef CGAL::Circulator_from_iterator<I>	Circulator;

typedef Ss::Face_const_iterator Face_const_iterator ;

template<class K>
void construct_total_graph( Polygon_2& poly, CGAL::Straight_skeleton_2<K> const& interior_ss, CGAL::Straight_skeleton_2<K> const& exterior_ss, std::vector<BridgingGraph>& bg)
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

		std::cout << "bg iss :" << bgNode.fid_iss << " bg ess: " << bgNode.fid_ess << std::endl;

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
	int same_v;
	double epsilon = 1e-16;

	Point ea(seg.vertex(0));
	Point eb(seg.vertex(1));

	for(Halfedge_const_iterator h = ss.halfedges_begin(); h != ss.halfedges_end(); ++h)
	{	
		same_v = 0;		//The number of same vertices

		//Find if the edge is identical to halfedge in straight skeleton structure
		Point ha(h->vertex()->point().x(), h->vertex()->point().y());
		Point hb(h->opposite()->vertex()->point().x(), h->opposite()->vertex()->point().y());

		if(abs((ha-ea).squared_length()) < epsilon )
			same_v ++;
		if(abs((ha-eb).squared_length()) < epsilon )
			same_v ++;
		if(abs((hb-ea).squared_length()) < epsilon )
			same_v ++;
		if(abs((hb-eb).squared_length()) < epsilon )
			same_v ++;

		if(same_v == 2)	//if they are identical
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
	//std::cout << "Face #" << f->id() << ": "<< f->halfedge()->vertex()->point().x() << " " <<f->halfedge()->vertex()->point().y() <<", "<< f->halfedge()->opposite()->vertex()->point().x() <<" "<<f->halfedge()->opposite()->vertex()->point().y() << std::endl;
}
template<class K>
void search_bridgegraph(std::vector<BridgingGraph> const& bg_list, Segment const& e, BridgingGraph& bg)
{
	//***** This function can be improved, because of special order in CGAL datastructure
	int same_v;
	double epsilon = 1e-16;

	Point ea(e.vertex(0));
	Point eb(e.vertex(1));

	for(I b=bg_list.begin(); b != bg_list.end(); ++b)
	{
		same_v = 0;
		Segment s = b->cut_edge;
		Point sa(s.vertex(0));
		Point sb(s.vertex(1));

		if(abs((ea-sa).squared_length()) < epsilon) ++same_v;
		if(abs((ea-sb).squared_length()) < epsilon) ++same_v;
		if(abs((eb-sa).squared_length()) < epsilon) ++same_v;
		if(abs((eb-sb).squared_length()) < epsilon) ++same_v;

		if(same_v == 2) { bg = b; return;}
	}
}

template<class K>
void search_bridgegraph(std::vector<BridgingGraph> const& bg_list, int const& fid, int const& offset, BridgingGraph& bg)
{
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
				return;
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
				return;
			}
		}
	}
}
#endif