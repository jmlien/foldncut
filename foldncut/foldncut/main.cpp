#include <boost/shared_ptr.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/create_straight_skeleton_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Straight_skeleton_2<K> Ss;
typedef boost::shared_ptr<Ss> SsPtr;


template<class K>
void print_point ( CGAL::Point_2<K> const& p )
{
  std::cout << "(" << p.x() << "," << p.y() << ")" ;
}

template<class K>
void print_straight_skeleton( CGAL::Straight_skeleton_2<K> const& ss )
{
  typedef CGAL::Straight_skeleton_2<K> Ss ;

  typedef typename Ss::Vertex_const_handle     Vertex_const_handle ;
  typedef typename Ss::Halfedge_const_handle   Halfedge_const_handle ;
  typedef typename Ss::Halfedge_const_iterator Halfedge_const_iterator ;

  Halfedge_const_handle null_halfedge ;
  Vertex_const_handle   null_vertex ;

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

int main()
{
	Polygon_2 poly;
	poly.push_back(Point(-1, -1));
	poly.push_back(Point(0, -12));
	poly.push_back(Point(1, -1));
	poly.push_back(Point(12, 0));
	poly.push_back(Point(1, 1));
	poly.push_back(Point(0, 12));
	poly.push_back(Point(-1, 1));
	poly.push_back(Point(-12, 0));

	//Via iterator pair
	SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());

	//Directy wat
	double IMaxOffset = 5;
	SsPtr oss = CGAL::create_exterior_straight_skeleton_2(IMaxOffset, poly);

	print_straight_skeleton(*iss);
	print_straight_skeleton(*oss);
}