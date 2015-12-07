//Straight Skeleton
#include <boost/shared_ptr.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/create_straight_skeleton_2.h>
#include <qline.h>

#include <CGAL/circulator.h>

//For Straight Skeleton
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Straight_skeleton_2<K> Ss;
typedef boost::shared_ptr<Ss> SsPtr;

//For Stitched Graph
typedef K::Segment_2 Segment;
typedef Polygon_2::Edge_const_circulator Edge_const_circulator;
typedef Ss::Halfedge_const_iterator Halfedge_const_iterator ;
typedef Ss::Halfedge_handle Halfedge_handle;
typedef struct stitchingGraph
{
	Segment cut_edge;		
	int fid_iss;
	int fid_ess;

}TOTALGRAPH;
typedef std::vector<stitchingGraph>::iterator I;
typedef CGAL::Circulator_from_iterator<I>	Circulator;

typedef Ss::Face_const_iterator Face_const_iterator ;

//For Visualization
#define MAX_DOUBLE (std::numeric_limits<double>::max)()
#define MIN_DOUBLE std::numeric_limits<double>::denorm_min()

template<class K>
void print_point ( CGAL::Point_2<K> const& p )
{
	std::cout << "(" << p.x() << "," << p.y() << ")" ;
}

//This function simply shows the bisectors and contours of straight skeleton in command window. 
template<class K>
void print_straight_skeleton( CGAL::Straight_skeleton_2<K> const& ss )
{

	typedef typename Ss::Halfedge_const_iterator Halfedge_const_iterator ;

	std::cout << "Straight skeleton with " << ss.size_of_vertices()
		<< " vertices, " << ss.size_of_halfedges()
		<< " halfedges and " << ss.size_of_faces()
		<< " faces" << std::endl ;

	for ( Halfedge_const_iterator i = ss.halfedges_begin(); i != ss.halfedges_end(); ++i )
	{
		print_point(i->opposite()->vertex()->point()) ;
		std::cout << "->" ;
		print_point(i->vertex()->point());
		std::cout << " " << ( i->is_bisector() ? "bisector" : "contour" ) << std::endl;
	}
}

//This function extracts the output from straight skeleton in CGAL and converts data to Qt type. 
template<class K>
void convert_straight_skeleton( CGAL::Straight_skeleton_2<K> const& ss , std::list<QLineF> &bisectors)
{	
	Point opp, p;

	for ( Halfedge_const_iterator i = ss.halfedges_begin(); i != ss.halfedges_end(); ++i )
	{
		if(i->is_bisector()){
			opp = i->opposite()->vertex()->point();
			p = i->vertex()->point();
			bisectors.push_back(QLineF(opp.x(), opp.y(), p.x(), p.y()));
		}
	} 
}
template<class K>
void construct_total_graph( Polygon_2& poly, CGAL::Straight_skeleton_2<K> const& interior_ss, CGAL::Straight_skeleton_2<K> const& exterior_ss, std::vector<stitchingGraph>& tg)
{
	if(poly.is_empty()) {std::cout << "Polygon is empty. Please check your model again."<< std::endl; exit(-1);}

	//Make circulator in the same order with polygon
	//This circulator contains information about cut edges and connected face to cut edges
	Edge_const_circulator root =  poly.edges_circulator();
	Edge_const_circulator e = poly.edges_circulator();
	
	do{
		stitchingGraph sgNode = stitchingGraph();
		sgNode.cut_edge = Segment(Point(e->vertex(0).x(), e->vertex(0).y()), Point(e->vertex(1).x(), e->vertex(1).y()));
		sgNode.fid_iss = find_corresponding_info(interior_ss, sgNode.cut_edge, 0);
		sgNode.fid_ess = find_corresponding_info(exterior_ss, sgNode.cut_edge, interior_ss.size_of_faces());

		if(sgNode.fid_iss == -1) {std::cout << "Cannot find match cutedges in straight skeleton structure" << std::endl; exit(-1);}
		if(sgNode.fid_ess == -1) {std::cout << "Cannot find match cutedges in straight skeleton structure" << std::endl; exit(-1);}
		tg.push_back(sgNode);
		++e;

	}while( e != root);
	
}
template<class K>
int find_corresponding_info(CGAL::Straight_skeleton_2<K> const& ss, Segment& seg, int offset)
{
	//***** This function can be improved, becuase both polygon and straight skeleton structures are stored in particular order. 
	int same_v;
	double epsilon = 1e-16;

	for(Halfedge_const_iterator h = ss.halfedges_begin(); h != ss.halfedges_end(); ++h)
	{	
		same_v = 0;		//The number of same vertices

		//Find if the edge is identical to halfedge in straight skeleton structure
		Point ha(h->vertex()->point().x(), h->vertex()->point().y());
		Point hb(h->opposite()->vertex()->point().x(), h->opposite()->vertex()->point().y());
		Point ea(seg.vertex(0).x(), seg.vertex(0).y());
		Point eb(seg.vertex(1).x(), seg.vertex(1).y());

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
			if(h->face() != NULL) return h->face()->id();
			else if(h->opposite()->face() != NULL) return h->opposite()->face()->id();
			else return -1;
		}
	}
	return -1;
}

//This function can help move on from cut edges to next skeleton faces(where reflected ray exists). 
template<class K>
void find_skeleton_face(CGAL::Straight_skeleton_2<K> const& ss, int fid, int offset)
{
	Face_const_iterator f = ss.faces_begin();
	int n = fid - offset;
	while( n != 0 )
	{
		++f;
		--n;
	}
	//std::cout << "Face #" << f->id() << ": "<< f->halfedge()->vertex()->point().x() << " " <<f->halfedge()->vertex()->point().y() <<", "<< f->halfedge()->opposite()->vertex()->point().x() <<" "<<f->halfedge()->opposite()->vertex()->point().y() << std::endl;
}