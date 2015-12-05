#include "handling_SS.h"

#include <iostream>
#include <fstream>

//Graphical View
#include <qapplication.h>
#include <qgraphicsscene.h>
#include <qgraphicsview.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>
#include <CGAL/number_utils.h>

#define MAX_DOUBLE (std::numeric_limits<double>::max)()
#define MIN_DOUBLE std::numeric_limits<double>::denorm_min()
#define PAPER_THRESHOLD 10.0

qreal minX = MAX_DOUBLE; qreal minY = MAX_DOUBLE;
qreal maxX = MIN_DOUBLE; qreal maxY = MIN_DOUBLE;
qreal width = 0;
qreal height = 0;

void read_file(char* filename, Polygon_2& poly, QPolygonF& qt_poly);
void createQTscene(QGraphicsScene& s, QPolygonF& qt_poly, std::list<QLineF>& qt_bis);

int main(int argc, char** argv)
{
	//Reading polygon files
	Polygon_2 poly; QPolygonF qt_polygon;	
	read_file("square.txt", poly, qt_polygon);

	//Using straight skeleton via iterator pair
	SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());

	//Extracts skeleton from cgal to qt
	std::list<QLineF> bis;
	convert_straight_skeleton(*iss, bis);

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

void read_file(char* filename, Polygon_2& poly, QPolygonF& qt_poly)
{
	std::ifstream f; 
	double x, y;

	f.open(filename);
	if(!f.is_open()) { std::cout << "Input polygon file is missing" << std::endl; exit(-1);}

	while(!f.eof())
	{
		f >> x >> y ;
		
		if(x < minX) minX = x; 		if(x > maxX) maxX = x;
		if(y < minY) minY = y;		if(y > maxY) maxY = y;

		poly.push_back(Point(x, y));
		qt_poly.push_back(QPointF(x, y));
	}

	width = maxX - minX;
	height = maxY - minY;

	if(width == 0 && height == 0) { std::cout << "Input polygon has no volume. Please check the file again" << std::endl; } 

	f.close();
}

void createQTscene(QGraphicsScene& s, QPolygonF& qt_poly, std::list<QLineF>& qt_bis)
{
	//Set the coordinates
	s.setSceneRect(minX-PAPER_THRESHOLD, minY-PAPER_THRESHOLD, width+2*PAPER_THRESHOLD, height+2*PAPER_THRESHOLD);	

	//Draw paper
	QPolygonF paper;															
	paper.push_back(QPointF(minX-PAPER_THRESHOLD, minY-PAPER_THRESHOLD));
	paper.push_back(QPointF(minX-PAPER_THRESHOLD, maxY+PAPER_THRESHOLD));
	paper.push_back(QPointF(maxX+PAPER_THRESHOLD, maxY+PAPER_THRESHOLD));
	paper.push_back(QPointF(maxX+PAPER_THRESHOLD, minY-PAPER_THRESHOLD));
	s.addPolygon(paper, QPen(QBrush(Qt::SolidLine), width/100.0));

	//Draw polygon
	s.addPolygon(qt_poly, QPen(QBrush(Qt::SolidLine), width/100.0));		

	//Draw skeleton
	QLineF line;
	while(qt_bis.size() !=0 )														
	{
		line = qt_bis.front();
		s.addLine(line, QPen(QBrush(Qt::blue), width/100.0));
		qt_bis.pop_front();
	}
}

//Using straight skeleton directly
//double IMaxOffset = 5;
//SsPtr oss = CGAL::create_exterior_straight_skeleton_2(IMaxOffset, poly);
