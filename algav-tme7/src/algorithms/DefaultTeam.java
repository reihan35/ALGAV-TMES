package algorithms;

import java.awt.Point;
import java.util.ArrayList;

/***************************************************************
 * TME 1: calcul de diamètre et de cercle couvrant minimum.    *
 *   - Trouver deux points les plus éloignés d'un ensemble de  *
 *     points donné en entrée.                                 *
 *   - Couvrir l'ensemble de poitns donné en entrée par un     *
 *     cercle de rayon minimum.                                *
 *                                                             *
 * class Circle:                                               *
 *   - Circle(Point c, int r) constructs a new circle          *
 *     centered at c with radius r.                            *
 *   - Point getCenter() returns the center point.             *
 *   - int getRadius() returns the circle radius.              *
 *                                                             *
 * class Line:                                                 *
 *   - Line(Point p, Point q) constructs a new line            *
 *     starting at p ending at q.                              *
 *   - Point getP() returns one of the two end points.         *
 *   - Point getQ() returns the other end point.               *
 ***************************************************************/
import supportGUI.Circle;
import supportGUI.Line;

public class DefaultTeam {

  // calculDiametre: ArrayList<Point> --> Line
  //   renvoie une paire de points de la liste, de distance maximum.
  public Line calculDiametre(ArrayList<Point> points) {
    if (points.size()<3) {
      return null;
    }

    Point p=points.get(0);
    Point q=points.get(1);
    double dist = p.distance(q);
    
    
    for (int i = 0; i < points.size(); i++) {
    	Point ptI = points.get(i);
    	for(int j = i; j < points.size(); j++) {
    		Point ptJ = points.get(j);
    		double newDist = ptI.distance(ptJ);
    		if (newDist > dist) {
    			p = ptI;
    			q = ptJ;
    			dist = newDist;
    		}
    	}
    }

    return new Line(p,q);
  }
  
  public Point calculDistMax(Point p, ArrayList<Point> points) {
	    if (points.size()<3) {
	      return null;
	    }

	    Point q=points.get(0);
	    double dist = p.distance(q);
	    
	    for (Point p1 : points) {
	    	double newDist = p1.distance(p);
	    	if (newDist > dist) {
	   			q = p1;
	   			dist = newDist;
	   		}
    	}
	    return q;
	  }
  
  public int[] media(Point p, Point q) {
	Point midpoint = new Point((p.x + q.x)/2,(p.y + q.y)/2);
	float slope = (q.y - p.y)/(q.x - p.x);
	float m = -1/slope;
	float b = midpoint.y - m*midpoint.x;
	
	int res[] = {(int)m, (int)b};
	return res;
	
  }

  public Circle cercleCirconscrit(Point p, Point q, Point r) {
	  int med1[] = media(p, q);
	  int med2[] = media(q, r);
	  int coeffX = med1[0] - med2[0];
	  int bSum = med2[1] - med1[0];
	  int x = bSum / coeffX;
	  int y = med1[0]*x + med1[0];
	  Point center = new Point(x, y);
	  
	  double radius = center.distance(p);
	  
	  return new Circle(center, (int) radius);
	  
	  
  }
  // calculCercleMin: ArrayList<Point> --> Circle
  //   renvoie un cercle couvrant tout point de la liste, de rayon minimum.
//  public Circle calculCercleMin(ArrayList<Point> points) {
//    if (points.isEmpty()) {
//      return null;
//    }
//    
//
//    Point center=points.get(0);
//    int radius=100;
//    Line res = calculDiametre(points);
//
//    center.x = (res.getP().x + res.getQ().x)/2;
//    center.y = (res.getP().y + res.getQ().y)/2;
//    radius = (int) (res.getP().distance(res.getQ())/2);
//    boolean iscmin = true;
//    
//    for (Point p1 : points) {
//    	if (p1.distance(center)>radius) {
//    		iscmin = false;
//    		break;
//    	}
//    }
//    
//
//	System.out.println("salut");
//    if (!iscmin) {
//        Circle c = new Circle(center, 1000000000);
//    	for (int i = 0; i < points.size(); i++) {
//    		Point p = points.get(i);
//    		for(int j = i+1; j < points.size(); j++) {
//        		Point q = points.get(j);
//    			for(int k = j+1; k < points.size(); k++) {
//    	    		Point r = points.get(k);
//    	    		Circle c2 = cercleCirconscrit(p, q, r);
//
//    	    		for (Point p1 : points) {
//    	    	    	if (p1.distance(center)>c2.getRadius()) {
//    	    	    		iscmin = false;
//    	    	    		break;
//    	    	    	}
//    	    	    }
//    	    		if (iscmin && c2.getRadius() < c.getRadius()) {
//    	    			c = c2;
//    	    		}
//
//    			}
//    		}
//    	}
//    }
//    return new Circle(center,radius);
//  }
  
  public int crossProduct(Point a, Point b) {
	  return a.x * b.y - a.y*b.x;
  }
  
  // p est dans le triangle ABC? avec produit vectoriel, non teste
//  public boolean isInTriangle(Point x, Point A, Point B, Point C) {
//	  
//	  if((crossProduct(B, x) <= 0 && crossProduct(B, C) <= 0) ||
//			 (crossProduct(B, x) >= 0 && crossProduct(B, C) >= 0)) {
//		  if((crossProduct(C, x) <= 0 && crossProduct(C, A) <= 0) ||
//				  (crossProduct(C, x) >= 0 && crossProduct(C, A) >= 0)) {
//			  if((crossProduct(C, x) <= 0 && crossProduct(C, B) <= 0) ||
//					  (crossProduct(C, x) >= 0 && crossProduct(C, B) >= 0)) {
//				  return true;
//			  }
//		  }
//	  }
//	  return false;
//  }
  
  public boolean isInTriangle(Point X, Point A, Point B, Point C) {
	  float l1 = ((float)((B.y - C.y) * (X.x - C.x) + (C.x - B.x) * (X.y - C.y)) )/((B.y - C.y) * (A.x - C.x) + (C.x - B.x) * (A.y - C.y));
	  float l2 = ((float)((C.y - A.y) * (X.x - C.x) + (A.x - C.x) * (X.y - C.y)))/((B.y - C.y) * (A.x - C.x) + (C.x - B.x) * (A.y - C.y));
	  double l3 = 1 - l1 - l2;
	  return l1>=0 && l1<=1 && l2>=0 && l2<=1 && l3>=0 && l3<=1;
  }
	  
  
  public ArrayList<Point> filtrage(ArrayList<Point> points){
	  ArrayList<Point> res = new ArrayList<>();
	  Point p = points.get(0);
	  int xmin, ymin, xmax, ymax;
	  
	  Point pN = p, pS = p, pO = p, pE = p;
	  xmin = p.x;
	  xmax = p.x;
	  ymin = p.y;
	  ymax = p.y;
	  
	  for (Point pt:points) {
		  int x = pt.x;
		  int y = pt.y;
		  if(x < xmin) {
			  pO = pt;
			  xmin = x;
		  }
		  if(x > xmax) {
			  pE = pt;
			  xmax = x;
		  }
		  if(y < ymin) {
			  pN = pt;
			  ymin = y;
		  }
		  if(y > ymax) {
			  pS = pt;
			  ymax = y;
		  }
	  }
	  
	  res.add(pN);
	  res.add(pE);
	  res.add(pS);
	  res.add(pO);
	  
	  for (Point pt:points) {
		  if (!(isInTriangle(pt, pN, pE, pO) || isInTriangle(pt, pS, pE, pO))){
			  res.add(pt);
		  }
	  }
	  
	  return res;
  }
  
//  public Point[] quickhull(ArrayList<Point> points) {
//	  Point p = points.get(0);
//	  int xmin, ymin, xmax, ymax;
//	  ArrayList<Point> res = new ArrayList<>(); 
//	  
//	  Point pN = p, pS = p, pO = p, pE = p;
//	  xmin = p.x;
//	  xmax = p.x;
//	  ymin = p.y;
//	  ymax = p.y;
//	  
//	  for (Point pt:points) {
//		  int x = pt.x;
//		  int y = pt.y;
//		  if(x < xmin) {
//			  pO = pt;
//			  xmin = x;
//		  }
//		  if(x > xmax) {
//			  pE = pt;
//			  xmax = x;
//		  }
//		  if(y < ymin) {
//			  pN = pt;
//			  ymin = y;
//		  }
//		  if(y > ymax) {
//			  pS = pt;
//			  ymax = y;
//		  }
//	  }
//	  
//	  res.add(pN);
//	  res.add(pE);
//	  res.add(pS);
//	  res.add(pO);
//	  
//	  for (Point pt:points) {
//		  
//	  }
//  }
	  
	  
  public Circle calculCercleMin(ArrayList<Point> points2) {
	  
	  
	  
	  //ArrayList<Point> points = (ArrayList<Point>) points2.clone();
	  ArrayList<Point> points = filtrage(points2);
	  Point dummy=points.get((int)Math.random()*points2.size());
	  Point p = calculDistMax(dummy, points);
	  Point q = calculDistMax(p, points);
	  Point c = new Point((p.x+q.x)/2, (p.y+q.y)/2);
	  Circle cercle = new Circle(c, (int) c.distance(p));
	  points.remove(p);
	  points.remove(q);
	  
	  while(!points.isEmpty()) {
		  Point s = points.get(0);
		  System.out.println(points.size());
		  if (s.distance(c) < cercle.getRadius()) {
			  points.remove(s);
		  }
		  else {
			  Point newC = new Point();
			  double newR = (c.distance(s) + cercle.getRadius())/ 2.0;
			  double alpha = newR/c.distance(s);
			  double beta = (c.distance(s)-newR)/c.distance(s);
			  newC.x = (int) (alpha * c.x + beta * s.x);
			  newC.y = (int) (alpha * c.y + beta * s.y);
			  c = newC;
			  cercle = new Circle(newC, (int) newR+1);
		  }
		  
		  
	  }
	  return cercle;
  }
}
