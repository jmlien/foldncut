/*
	Project : fold n cut

	This project works with polygon without holes. 


	Dependency : CGAL 4.7 (libQGLViewer 2.6.3)
				 Qt 5.5 
*/

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
	construct_bridging_graph(poly, *iss, *ess, bg);

	//Compute perpendiculars
	std::list<Perpendiculars> ppd;
	generate_perpendiculars(*iss, *ess, bg, ppd);
	//deduplicate_perpendiculars<K>(ppd);

	//Assign mountain and valley
	std::list<Segment> mt;
	std::list<Segment> vl;
	MountainValley<K>(*iss, *ess, bg, ppd,  mt, vl);
	//unique_perpendiculars(ppd);

	//Extracts skeleton from cgal and converts them to qt
	std::list<QLineF> bis_qt;
	convert_straight_skeleton(*iss, bis_qt);
	convert_straight_skeleton(*ess, bis_qt);

	std::list<QLineF> ppd_qt;
	convert_perpendiculars<K>(ppd, ppd_qt);
	
	std::list<QLineF> mt_qt, vl_qt;
	convert_mountain_valley<K>(mt, mt_qt);
	convert_mountain_valley<K>(vl, vl_qt);
	
	//SETUP QT 
	//Create applicaiton
	QApplication app(argc, argv);

	//Create scene
	QGraphicsScene scene;
	createQTscene(scene, poly_qt, bis_qt, ppd_qt, mt_qt, vl_qt);

	//Create view
	QGraphicsView* view = new QGraphicsView(&scene);
	CGAL::Qt::GraphicsViewNavigation nav;
	view->installEventFilter(&nav);
	view->viewport()->installEventFilter(&nav);
	view->setRenderHint(QPainter::Antialiasing);
	//view->keyPressEvent();
	
	//Show view
	view->show();

	return app.exec();
}