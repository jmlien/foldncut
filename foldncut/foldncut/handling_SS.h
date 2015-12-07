//Straight Skeleton
#include <boost/shared_ptr.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/create_straight_skeleton_2.h>
#include <qline.h>

//CGAL : Straight Skeleton
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Straight_skeleton_2<K> Ss;
typedef boost::shared_ptr<Ss> SsPtr;

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

	//typedef typename Ss::Vertex_const_handle     Vertex_const_handle ;
	//typedef typename Ss::Halfedge_const_handle   Halfedge_const_handle ;
	typedef typename Ss::Halfedge_const_iterator Halfedge_const_iterator ;

	//Halfedge_const_handle null_halfedge ;
	//Vertex_const_handle   null_vertex ;

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
	//typedef typename Ss::Vertex_const_handle     Vertex_const_handle ;
	//typedef typename Ss::Halfedge_const_handle   Halfedge_const_handle ;
	typedef typename Ss::Halfedge_const_iterator Halfedge_const_iterator ;

	//Halfedge_const_handle null_halfedge ;
	//Vertex_const_handle   null_vertex ;
	CGAL::Point_2<K> opp, p;

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
void construct_total_graph( CGAL::Straight_skeleton_2<K> const& interior_ss, CGAL::Straight_skeleton_2<K> const& exterior_ss)
{
	//*****Try to use triedge, but doens't contain information of edges
	//*****Try to use circulator, failed
	typedef typename Ss::Vertex_const_iterator	Vertex_const_iterator;
	typedef typename Ss::Halfedge_const_handle   Halfedge_const_handle ;
	typedef typename Ss::Halfedge_const_iterator Halfedge_const_iterator ;

	Halfedge_const_iterator h;
	//Halfedge_around_vertex_const_circulator e;
	for ( Vertex_const_iterator v = interior_ss.vertices_begin(); v != interior_ss.vertices_end(); ++v )
	{
		std::cout <<"vertex:"<< v->point().x() <<" "<< v->point().y() << std::endl; 
		h = v->halfedge_around_vertex_begin();
		if(v->is_contour()) std::cout << "on contour" << std::endl;
		if(v->is_skeleton()) std::cout << "skeleton" << std::endl;
		
		std::cout << std::endl;
	}
}

void find_corresponding_cutedge(/* skeleton vertex*/ /* cut edges*/)
{

}