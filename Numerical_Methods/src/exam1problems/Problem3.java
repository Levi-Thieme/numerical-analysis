package exam1problems;

import java.util.function.Function;

import numerical_functions.NumericalFunction;

public class Problem3 {

	public static void main(String[] args) {
		//[0, 1.25], x = {.25, .50, 1.0}, y = {25.2, 49.2, 96.4}
		Double[] x = {.25, .5, 1.0};
		Double[] y = {25.2, 49.2, 96.4};
		double[][] dividedDiffs = NumericalFunction.dividedDiff(x, y);
		printDifferencesTriangle(dividedDiffs);
		
		Function<Double, Double> f = NumericalFunction.newtonInterp(x, y);
		System.out.println(f.apply(.75));
	}
	
	/**
	 * Helper function for printing the divided differences matrix.
	 * 
	 * @param dividedDiffs A double matrix containing divided differences.
	 */
	private static void printDifferencesTriangle(double[][] dividedDiffs) {
		for (double[] row : dividedDiffs) {
			for (double rowColumn : row) {
				if (rowColumn != 0.0) {
					System.out.format("%6.2f  ", rowColumn);
				}
			}
			System.out.println("\n");
		}
	}

}
