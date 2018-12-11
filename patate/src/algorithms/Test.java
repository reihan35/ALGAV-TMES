package algorithms;

import java.awt.Point;

public class Test {

	public static void main(String[] args) {
		DefaultTeam de = new DefaultTeam();
		Point A = new Point(586, 259);
		Point B = new Point(401, 442);
		Point C = new Point(1, 1);
		Point D = new Point(3, 2);	
		Tuple<Double, Double> t= de.droite(C,  D);
		System.out.println(de.droite_dist(A, A, B));
		System.out.println(t.getX() + " " + t.getY());

	}

}
