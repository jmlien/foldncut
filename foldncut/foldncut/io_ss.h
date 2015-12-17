#ifndef SS_IO_H
#define SS_IO_H
#include "globals.h"

#include <iostream>
#include <fstream>

#include <qline.h>

//Graphical View
#include <qapplication.h>
#include <qgraphicsscene.h>
#include <qgraphicsview.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>
#include <CGAL/number_utils.h>

qreal width = 0;
qreal height = 0;

void read_file(char* filename, Polygon_2& poly, QPolygonF& qt_poly)
{
	//*****y coordinate is flipped!
	std::ifstream f; 
	double x, y;

	f.open(filename);
	if(!f.is_open()) { std::cout << "Input polygon file is missing" << std::endl; exit(-1);}
	else std::cout << "Opened file successfully" << std::endl;

	minX = MAX_DOUBLE; minY = MAX_DOUBLE;
	maxX = MIN_DOUBLE; maxY = MIN_DOUBLE;

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

void createQTscene(QGraphicsScene& s, QPolygonF& poly, std::list<QLineF>& bis, std::list<QLineF>& ppd, std::list<QLineF>& mt, std::list<QLineF>& vl)
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
	s.addPolygon(poly, QPen(QBrush(Qt::SolidLine), width/100.0));		

	//Draw skeleton
	QLineF line;
	while(bis.size() !=0 )														
	{
		line = bis.front();
		s.addLine(line, QPen(QBrush(Qt::red), width/300.0));
		bis.pop_front();
	}

	//Draw perpendiculars
	while(ppd.size() !=0 )														
	{
		line = ppd.front();
		s.addLine(line, QPen(QBrush(Qt::blue), width/300.0));
		ppd.pop_front();
	}

	//Draw mountain
	while(mt.size() !=0 )														
	{
		line = mt.front();
		s.addLine(line, QPen(QBrush(Qt::green), width/300.0));
		mt.pop_front();
	}

	//Draw valley
	while(vl.size() !=0 )														
	{
		line = vl.front();
		s.addLine(line, QPen(QBrush(Qt::yellow), width/300.0));
		vl.pop_front();
	}
}

template<class K>
void print_point ( CGAL::Point_2<K> const& p )
{
	std::cout << "(" << p.x() << "," << p.y() << ")" ;
}

//This function simply shows the bisectors and contours of straight skeleton in command window. 
template<class K>
void print_straight_skeleton_vertex( CGAL::Straight_skeleton_2<K> const& ss )
{

	typedef typename Ss::Vertex_const_iterator Vertex_const_iterator ;

	for ( Vertex_const_iterator i = ss.vertices_begin(); i != ss.vertices_end(); ++i )
	{
		std::cout << (i->is_contour() ? "contour" : "skeleton") <<std::endl;
		print_point(i->point());
		std::cout << std::endl;
	}
}

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
			CGAL::exact(opp.x());
			CGAL::exact(opp.y());
			CGAL::exact(p.x());
			CGAL::exact(p.y());
			bisectors.push_back(QLineF(opp.x(), opp.y(), p.x(), p.y()));
		}
	} 
}

template<class K>
void convert_perpendiculars(std::list<Perpendiculars>& ppd, std::list<QLineF> &perpendiculars)
{
	if(ppd.empty()) {std::cout << "There is no perpendiculars! Please check perpendiculars again."<< std::endl; return;}

	Point opp, p;

	Segment s;

	while(!ppd.empty())
	{
		s = ppd.front().seg;
		p = s.vertex(0);
		opp = s.vertex(1);
		perpendiculars.push_back(QLineF(opp.x(), opp.y(), p.x(), p.y()));
		ppd.pop_front();
	}
}

template<class K>
void convert_mountain_valley(std::list<Segment>& mtvl, std::list<QLineF> &mtvl_qt)
{
	if(mtvl.empty()) {std::cout << "There is no perpendiculars! Please check perpendiculars again."<< std::endl; return;}

	Point opp, p;

	Segment s;

	while(!mtvl.empty())
	{
		s = mtvl.front();
		p = s.vertex(0);
		opp = s.vertex(1);

		mtvl_qt.push_back(QLineF(opp.x(), opp.y(), p.x(), p.y()));

		mtvl.pop_front();
	}
}
#endif