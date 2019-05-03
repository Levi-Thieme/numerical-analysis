package linear_algebra;

import java.util.ArrayList;

public class Norm {
	
	public static double euclideanNorm(ArrayList<Double> v) {
		double norm = 0.0;
		for (double entry : v) {
			norm += entry * entry;
		}
		return Math.sqrt(norm);
	}
	
	public static double euclideanNorm(double[] v) {
		double norm = 0.0;
		for (double entry : v) {
			norm += entry * entry;
		}
		return Math.sqrt(norm);
	}
	
	public static double infiniteNorm(ArrayList<Double> v) {
		double norm = Math.abs(v.get(0));
		for (double entry : v) {
			norm = Math.max(norm, Math.abs(entry));
		}
		return norm;
	}
	public static double infiniteNorm(double[]  v) {
		double norm = Math.abs(v[0]);
		for (double entry : v) {
			norm = Math.max(norm, Math.abs(entry));
		}
		return norm;
	}
	
	/**
	 * Cauchy-Bunyakovsky-Schwarz Inequality for Sums
	 * 
	 * 
	 * @param x
	 * @param y
	 * @return
	 */
	public static double CBSInequalitySums(ArrayList<Double> x, ArrayList<Double> y) {
		return euclideanNorm(x) * euclideanNorm(y);
	}
	
	public static double distance(ArrayList<Double> x, ArrayList<Double> y) {
		double distance = 0.0;
		for (int i = 0; i < x.size(); i++) {
			double difference = x.get(i) - y.get(i);
			distance += difference * difference;
		}
		return Math.sqrt(distance);
	}
	
	public static double determinant2D(double[][] m) {
		//ad - bc
		return m[0][0] * m[1][1] - m[1][0] * m[0][1];
	}
	
	public static double matrixInfiniteNorm(double[][] m) {
		double norm = 0.0;
		double rowSum = 0.0;
		for (int i = 0; i < m.length; i++) {
			rowSum = 0.0;
			for (int j = 0; j < m[i].length; j++) {
				rowSum += Math.abs(m[i][j]);
			}
			norm = Math.max(norm, rowSum);
		}
		return norm;
	}
}


