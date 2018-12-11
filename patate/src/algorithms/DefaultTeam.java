package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;

import supportGUI.Circle;
import supportGUI.Line;

public class DefaultTeam {

	
	public Point c(Line d, Point q) {
		Tuple<Double, Double> ab = droite(d.getP(), d.getQ());
		//on calcule lequation de la parallele passant par q
		System.out.println(d.getP() + " " + d.getQ());
		System.out.println(q);
		//System.out.println(ab.getX() + " " + ab.getY());
		double a = ab.getX();
		double b = q.getY() - q.getX()*a;
		System.out.println("a vaut " + a + " b vaut " + b);
		//on calcule lequation de la perpendiculaire passant par A
		double aP = -1/ab.getX();
		double bP = d.getP().getY() - d.getP().getX()*aP;
		
		System.out.println("a vaut " + aP + " b vaut " + bP);
		
		double x = (bP-b)/(a-aP);
		double y = a*x + b;
		System.out.println("x : " + x + "y ! "+y);
		return new Point((int)x, (int)y);
		
	}
	
	public Point d(Line d, Point q) {
		Tuple<Double, Double> ab = droite(d.getP(), d.getQ());
		//on calcule lequation de la parallele passant par q

		double a = ab.getX();
		double b = q.getY() - q.getX()*a;
		//on calcule lequation de la perpendiculaire passant par A
		double aP = -1/ab.getX();
		double bP = d.getQ().getY() - d.getQ().getX()*aP;
		
		double x = (bP-b)/(a-aP);
		double y = a*x + b;
		return new Point((int)x, (int)y);
		
	}
	
	public ArrayList<Point> calculRectangleMin(ArrayList<Point> ev, Line diam) {
		Point A = diam.getP();
		Point B = diam.getQ();
		
		Point right = null;
		double maxDistR = 0;
		Point left = null;
		double maxDistL = 0;
		
		for(Point q: ev) {
			if(crossProduct(A, q, A, B) > 0) {
				double newDistR = droite_dist(q, A, B);
				if(newDistR > maxDistR) {
					maxDistR = newDistR;
					right = q;
				}
			}
			else {
				if(crossProduct(A, q, A, B) < 0) {
					double newDistL = droite_dist(q, A, B);
					if(newDistL > maxDistL) {
						maxDistL = newDistL;
						left = q;
					}
				}
			}
		}
		ArrayList<Point> res = new ArrayList<>();
		res.add(c(diam, right));
		res.add(d(diam, right));
		res.add(d(diam, left));
		res.add(c(diam, left));
		return res;
	}
	
	public ArrayList<Point> enveloppeConvexe(ArrayList<Point> points){
		return calculRectangleMin(enveloppeConvexeb(points), calculDiametre(points));
	}
	
	
  // calculDiametre: ArrayList<Point> --> Line
  //   renvoie une pair de points de la liste, de distance maximum.
  public Line calculDiametre(ArrayList<Point> points) {
	ArrayList<Point> ev = enveloppeConvexeb(points);
    Tuple<Point,Point> t = shamos(paire_antipodales(ev));
    return new Line(t.getX(),t.getY());
  }

  // calculDiametreOptimise: ArrayList<Point> --> Line
  //   renvoie une pair de points de la liste, de distance maximum.
  public Line calculDiametreOptimise(ArrayList<Point> points) {
	  
		ArrayList<Point> ev = enveloppeConvexe(points);
	    Tuple<Point,Point> t = shamos(paire_antipodales(ev));
	    //return calculRectangleMin(ev, new Line(t.getX(),t.getY()));
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


  
  
  public int crossProduct(Point p, Point q, Point s, Point t) {
		return ((q.x-p.x)*(t.y-s.y))-((q.y-p.y)*(t.x-s.x));
	}
  
  /*
  public int crossProduct(Point a, Point b) {
		return a.x * b.y - a.y*b.x;
	}*/
// enveloppeConvexe: ArrayList<Point> --> ArrayList<Point>
//   renvoie l'enveloppe convexe de la liste.
  /*
	public ArrayList<Point> enveloppeConvexe(ArrayList<Point> points){
	  if (points.size()<3) {
	    return null;
	  }
	  //points = filtre_pix(points); 
	  points = filtrage(points);
	  ArrayList<Point> enveloppe = new ArrayList<Point>();
	
	  boolean isCote = false;
	  
	  for (Point A : points) {
	
		  for(Point B : points) {
	  		if(B.equals(A))
	  			continue;
	  		int prodVecPrec = 0;
	  		Point R = B;
	  		for (Point rr:points) {
	  			if(crossProduct(A, B, A, rr) != 0) {
	  				R = rr;
	  				break;
	  			}
	  		}
	  		double signeR = crossProduct(A, B, A, R);
	  		isCote = true;
	  		for(Point C: points) {
	  			double prodVec = crossProduct(A, B, A, C);
	  			System.out.println(prodVec*prodVecPrec);
	  			if (signeR * prodVec < 0) {
	  				isCote = false;
	  				break;
	  			}
	  			
	  		}
	  		if(isCote) {
	
	  			enveloppe.add(A);
	  			enveloppe.add(B);
	  		}
	
	  	}
	  }
	  
	  	
	  return enveloppe;
	}
	*/
	
	public ArrayList<Point> filtre_pix( ArrayList<Point> points){
		
		int cpt = 0;
		ArrayList<Point> points2 = new ArrayList<>();
		
		
		Point min;
		Point max;
		
		for (Point p : points) {
			min = p;
			max = p;
			for (Point p2 : points ) {
				if(p.getX() == p2.getX()) {
					if (p2.getY() < min.getY()) {
						min = p;
					}
					
					if (p2.getY() > max.getY()) {
						max = p;
					}
				}
			}
			if (min.equals(max)){
				points2.add(min);
			}
			
			else {
				points2.add(min);
				points2.add(max);
			}
		}
		return points2;
	}
	  
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
	
	public Point trouverCote(Point p, ArrayList<Point> points) {
		boolean isCote;
		for(Point q : points) {
			if(q.equals(p))
	  			continue;
	  		Point r = p;
	  		for (Point rr:points) {
	  			if(crossProduct(p, q, p, rr) != 0) {
	  				r = rr;
	  				break;
	  			}
	  		}
	  		double signeR = crossProduct(p, q, p, r);
	  		isCote = true;
	  		for(Point s: points) {
	  			double prodVec = crossProduct(p, q, p, s);
	  			if (signeR * prodVec < 0) {
	  				isCote = false;
	  				break;
	  			}
	  		}
	  		if(isCote) {
	  			return q;
	  		}
		}
		throw new RuntimeException("Le point donne en parametre n'appartient pas a l'enveloppe convexe");
	}
	
	public double calculAngle(Point a,Point b,Point c, Point d) {
		double ux = b.getX() - a.getX(); 
		double uy = b.getY() - a.getY();
		double vx = d.getX() - c.getX();
		double vy = d.getY() - c.getY();
		
		double uv = (ux*vx + uy*vy);
		double normeu = Math.sqrt(ux*ux + uy*uy);
		double normev = Math.sqrt(vx*vx + vy*vy);
		
		return Math.acos(uv/(normeu*normev));
		
	}
	
	public Point pointAngleMin(Point p, Point q, ArrayList<Point> points) {
		double angleMin = 360;
		Point res = null;
		for(Point r:points) {
			double angle = calculAngle(p, q, q, r);
			System.out.println("l'angle vaut " + angle);
			if(angle < angleMin) {
				angleMin = angle;
				res = r;
			}
		}
		return res;
	}
	//JARVIS
//	public ArrayList<Point> enveloppeConvexe(ArrayList<Point> points){
//		double min = points.get(0).x;
//		ArrayList<Point> enveloppe = new ArrayList<Point>();
//		Point minp = null;
//		for (Point p: points) {
//			if (p.getX() < min) {
//				min = p.getX();
//				minp = p;
//			}
//		}
//		Point p = minp;
//		Point deb = p;
//		System.out.println("sa");
//		Point q = trouverCote(p, points);
//		Point r = pointAngleMin(p, q, points);
//		enveloppe.add(p);
//		enveloppe.add(q);
//		enveloppe.add(r);
//		
//		while(!r.equals(deb)) {
//			p = q;
//			q = r;
//			r = pointAngleMin(p, q, points);
//			enveloppe.add(r);
//		}
//		
//		return enveloppe;
//	}
//	
	public ArrayList<Point> enveloppeConvexeb(ArrayList<Point> points){
		ArrayList<Point> enveloppe = new ArrayList<Point>();
		Point[] bucketMin;
		Point[] bucketMax;
		
		double xMin = points.get(0).getX();
		double xMax = points.get(0).getX();
		for (Point p: points) {
			if (p.getX() < xMin) {
				xMin = p.getX();
			}
			if (p.getX() > xMax) {
				xMax = p.getX();
			}
		}
		
		
		
		int width = (int) (xMax-xMin) + 1;
		bucketMin = new Point[width];
		bucketMax = new Point[width];
		
		for(int i = 0; i < bucketMin.length; i++) {
			bucketMin[i] = null;
			bucketMax[i] = null;
		}

		for (Point p: points) {
			
			int i = (int) (p.getX() - xMin);
			if(bucketMin[i] != null) {
				if(p.getY() < bucketMin[i].getY()) {
					bucketMin[i] = p;
				}
				if(p.getY()>bucketMax[i].getY()) {
					bucketMax[i] = p;
				}
			}
			else {
				bucketMin[i] = p;
				bucketMax[i] = p;
			}
		}

		ArrayList<Point> minMax = new ArrayList<>();
		for(int i = 0; i < bucketMin.length; i++) {
			if(bucketMin[i] !=null)
				minMax.add(bucketMin[i]);
		}
		for(int i = bucketMax.length-1; i >=0; i--) {
			if(bucketMax[i] != null)
				minMax.add(bucketMax[i]);
		}
		minMax.add(bucketMin[0]);
		
		int i = 0;
		int cpt = 0;
		
		//return minMax;
		
		int l = minMax.size();
		while (cpt < 2*l) {
			cpt++;
			Point A = minMax.get(i%minMax.size());
			Point B = minMax.get((i+1)%minMax.size());
			Point C = minMax.get((i+2)%minMax.size());

			if(crossProduct(A, B, B, C) >= 0) {
				i++;
			}
			else {
				minMax.remove((i+1)%minMax.size());
				i = (i-1 + minMax.size())%minMax.size();
			}
		}
		
		return minMax;
		
	}
	
	public Tuple<Double,Double> droite(Point p1, Point p2){
		System.out.println("p1  : " + p1 + " p 2 : " + p2);
		 if (p1 == p2) {
			 return new Tuple<>(0.0, 0.0);
		 }
		 double a = ((double)(p2.y - p1.y))/(p2.x - p1.x);
		 double b = p1.y - (a * p1.x);
		 
		 return new Tuple<>(a,b);
	}
	
	public double droite_dist(Point p1, Point p2, Point p3) {
		Tuple<Double, Double> d = droite(p2,p3);
		if(d.getX()==0 && d.getY()==0)return 0;
		Double a = d.getX();
		Double b = d.getY();
		double dist = Math.abs(p1.y - a*p1.x - b) / (Math.sqrt(1 + b*b));
		return dist;
	}
	
	public HashSet<Tuple<Point,Point>> paire_antipodales(ArrayList<Point> points2){
		HashSet<Tuple<Point,Point>> A = new HashSet<>();
		System.out.println(points2);
		ArrayList<Point>points = new ArrayList<Point>(new LinkedHashSet<Point>(points2));
		System.out.println(points);
		int k = 1 ;
		int i;
		int j;
		int n = points.size();
		
		while(droite_dist(points.get(k),points.get(n-1),points.get(0)) < droite_dist(points.get(k+1),points.get(n-1),points.get(0))) {
			k = k + 1;
		}
		i = 0;
		j = k;
		System.out.println("k vaut " + k);

		while (i < k+1 && j < n) {
			Point pj = points.get(j%n);
			Point pi = points.get(i%n);
			Point pisuiv = points.get((i+1)%n);
			Point pjsuiv = points.get((j+1)%n);
			System.out.println("pi = " + pi + "pisuiv = " + pisuiv + "A = " + pj + " B = " + pjsuiv);
			while (droite_dist(pj,pi,pisuiv) < droite_dist(pjsuiv,pi,pisuiv) && j<n) {
				//System.out.println("sa1");
				//System.out.println("A = " + pj + " B = " + pjsuiv);
				A.add(new Tuple<Point,Point>(pi,pj));
				j++;
				pj = pjsuiv;
				pjsuiv = points.get((j+1)%n);
			}
			//System.out.println("sa");
			A.add(new Tuple<Point,Point>(points.get(i%n),points.get(j%n)));
			i++;
		}

		return A;
	}
	
	public Tuple<Point,Point> shamos(HashSet<Tuple<Point,Point>> paires_anti){
		
		double max = 0;
		Tuple<Point,Point> m = null; 
		

		for (Iterator<Tuple<Point, Point>> iterator = paires_anti.iterator(); iterator.hasNext();) {
			Tuple<Point, Point> p = iterator.next();
			System.out.println("p : "+p.getX()+p.getY());
			if (p.getX().distance(p.getY()) > max ) {
				max = p.getX().distance(p.getY());
				m = p;
			}
		}
		return m;
	}

}
