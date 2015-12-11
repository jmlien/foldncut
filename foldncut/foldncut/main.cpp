//#include "datastructure_ss.h"
#include "io_ss.h"
#include "foldncut.h"

int main(int argc, char** argv)
{
	//Reading polygon files
	Polygon_2 poly; QPolygonF qt_polygon;	
	read_file("models/bigstar.txt", poly, qt_polygon);

	//Compute straight skeleton
	//*****This program now works for polygon WITHOUT HOLES. 
	SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());
	SsPtr ess = CGAL::create_exterior_straight_skeleton_2(EMaxOffset, poly);
	
	//Constructing the graph combining the straight skeleton and polygon
	std::vector<stitchingGraph> stg;
	construct_total_graph(poly, *iss, *ess, stg);

	//Compute perpendiculars
	std::list<Halfedge_handle> ppd;
	generate_perpendiculars(*iss, *ess, stg, ppd);

	//Extracts skeleton from cgal and converts them to qt
	std::list<QLineF> bis;
	convert_straight_skeleton(*iss, bis);
	convert_straight_skeleton(*ess, bis);

	//SETUP QT 
	//Create applicaiton
	QApplication app(argc, argv);

	//Create scene
	QGraphicsScene scene;
	createQTscene(scene, qt_polygon, bis);

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