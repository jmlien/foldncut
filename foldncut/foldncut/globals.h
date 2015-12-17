#ifndef SS_GLOBALS_H
#define SS_GLOBALS_H

//Straight Skeleton
#include <boost/shared_ptr.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/create_straight_skeleton_2.h>

//For Straight Skeleton
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Line_2 Line;
typedef K::Segment_2 Segment;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Straight_skeleton_2<K> Ss;
typedef boost::shared_ptr<Ss> SsPtr;

//For paper (termainate condition)
#define MAX_DOUBLE (std::numeric_limits<double>::max)()
#define MIN_DOUBLE std::numeric_limits<double>::denorm_min()
#define PAPER_THRESHOLD 100.0

static double minX, minY, maxX, maxY;

typedef struct Perpendiculars
{
	Segment seg;
	int level;
	bool mountain;

}PERPENDICULARS;

/*

  foldncut.h  <-- datastructure_ss.h <-- globals.h
  
  ios_ss.h <-- globals.h

*/
#endif
