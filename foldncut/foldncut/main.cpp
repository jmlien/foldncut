
//Straight Skeleton
#include <boost/shared_ptr.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/create_straight_skeleton_2.h>

//Graphical View
#include <qapplication.h>
#include <qgraphicsscene.h>
#include <qgraphicsview.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>
//#include <CGAL/Qt/Converter.h>
#include <CGAL/number_utils.h>

//CGAL : Straight Skeleton
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Straight_skeleton_2<K> Ss;
typedef boost::shared_ptr<Ss> SsPtr;

//typedef CGAL::Qt::Converter<K> Converter; //Different from API on internet

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

template<class K>
void convert_straight_skeleton( CGAL::Straight_skeleton_2<K> const& ss , std::list<QLineF> &bisectors)
{
	typedef CGAL::Straight_skeleton_2<K> Ss ;

  typedef typename Ss::Vertex_const_handle     Vertex_const_handle ;
  typedef typename Ss::Halfedge_const_handle   Halfedge_const_handle ;
  typedef typename Ss::Halfedge_const_iterator Halfedge_const_iterator ;

  Halfedge_const_handle null_halfedge ;
  Vertex_const_handle   null_vertex ;
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

int main(int argc, char** argv)
{
	Polygon_2 poly; QPolygonF qt_polygon;	//Converter doesn't have Polygon Conversion from cgal to qt
	poly.push_back(Point(-1, -1)); qt_polygon.push_back(QPointF(-1, -1));
	poly.push_back(Point(0, -12)); qt_polygon.push_back(QPointF(0, -12));
	poly.push_back(Point(1, -1)); qt_polygon.push_back(QPointF(1, -1));
	poly.push_back(Point(12, 0)); qt_polygon.push_back(QPointF(12, 0));
	poly.push_back(Point(1, 1)); qt_polygon.push_back(QPointF(1, 1));
	poly.push_back(Point(0, 12)); qt_polygon.push_back(QPointF(0, 12));
	poly.push_back(Point(-1, 1)); qt_polygon.push_back(QPointF(-1, 1));
	poly.push_back(Point(-12, 0)); qt_polygon.push_back(QPointF(-12, 0));

	//Using straight skeleton via iterator pair
	SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());

	//Using straight skeleton directly
	double IMaxOffset = 5;
	SsPtr oss = CGAL::create_exterior_straight_skeleton_2(IMaxOffset, poly);

	//print_straight_skeleton(*iss);
	//print_straight_skeleton(*oss);

	//Convert skeletron from cgal to qt
	std::list<QLineF> bis;
	convert_straight_skeleton(*iss, bis);
	
	//SETUP QT 
	//Create applicaiton
	QApplication app(argc, argv);

	//Create scene
	QGraphicsScene scene;
	scene.setSceneRect(-15, -15, 30, 30);
	scene.addPolygon(qt_polygon, QPen(QBrush(Qt::SolidLine), 0.1)); 
	
	QLineF line;
	while(bis.size() !=0 )
	{
		line = bis.front();
		scene.addLine(line, QPen(QBrush(Qt::blue), 0.1));
		bis.pop_front();
	}
	
	//Create view
	QGraphicsView* view = new QGraphicsView(&scene);
	CGAL::Qt::GraphicsViewNavigation nav;
	view->installEventFilter(&nav);
	view->viewport()->installEventFilter(&nav);
	view->setRenderHint(QPainter::Antialiasing);

	//Show view
	view->show();

	return app.exec();
}