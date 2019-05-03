package exam1problems;

import numerical_functions.NumericalFunction;

public class Problem4 {

	public static void main(String[] args) {
		Double[] x = { 0.0, .25, .75, 1.0 };
		double[] y = { 1.0, Math.pow(Math.E, -.25), Math.pow(Math.E, -.75), 1.0 / Math.E };
		
		double[][] coefficients = NumericalFunction.clampedCubicSplineCoefficients(x, y, -1, -1/Math.E);
		printCubicSplines(x, coefficients);
	}
	
	private static void printCubicSplines(Double[] x, double[][] coefficients) {
		for (int i = 0; i < coefficients.length; i++) {
			double a = coefficients[0][i];
			double b = coefficients[1][i];
			double c = coefficients[2][i];
			double d = coefficients[3][i];
			//System.out.printf("(%8f, %8f, %8f, %8f)\n", a, b, c, d);
			System.out.printf("%8f + %8f(x - %2.2f) + %8f(x - %2.2f)^2 + %8f(x - %2.2f)^3\n", a, b, x[i], c, x[i], d, x[i]);
		}
	}

}
