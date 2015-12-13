#include "io_ss.h"
#include "creasepattern.h"

double EMaxOffset = PAPER_THRESHOLD * 1.5;

int main(int argc, char** argv)
{
	//Reading polygon files
	Polygon_2 poly; QPolygonF poly_qt;	
	read_file("models/turtle.txt", poly, poly_qt);

	//Compute straight skeleton
	//*****This program now works for polygon WITHOUT HOLES. 
	SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());
	SsPtr ess = CGAL::create_exterior_straight_skeleton_2(EMaxOffset, poly);
	
	//Constructing the graph combining the straight skeleton and polygon
	std::vector<BridgingGraph> bg;
	construct_total_graph(poly, *iss, *ess, bg);

	//Compute perpendiculars
	std::list<Segment> ppd;
	generate_perpendiculars(*iss, *ess, bg, ppd);

	//Extracts skeleton from cgal and converts them to qt
	std::list<QLineF> bis_qt;
	convert_straight_skeleton(*iss, bis_qt);
	convert_straight_skeleton(*ess, bis_qt);

	std::list<QLineF> ppd_qt;
	convert_perpendiculars<K>(ppd, ppd_qt);

	//SETUP QT 
	//Create applicaiton
	QApplication app(argc, argv);

	//Create scene
	QGraphicsScene scene;
	createQTscene(scene, poly_qt, bis_qt, ppd_qt);

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