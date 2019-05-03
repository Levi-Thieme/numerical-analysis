package tests;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.function.Function;

import org.junit.Test;

import numerical_functions.NumericalFunction;

public class CubicInterpTest {
	/**
	 * Tests a cubic spline approximation for f(x) = e^x
	 * 
	 */
	@Test
	public void test() {
		Function<Double, Double> f = (x) -> { return Math.pow(Math.E, x); };
		Double[] x = { 0.0, 1.0, 2.0, 3.0 };
		double[] y = { f.apply(x[0]), f.apply(x[1]), f.apply(x[2]), f.apply(x[3]) };
		
		double[][] coefficients = NumericalFunction.naturalCubicSplineCoefficients(x, y);
		ArrayList<Function<Double, Double>> S = NumericalFunction.cubicSplineInterp(x, coefficients, x.length);
		
		for (int i = 0; i < coefficients.length; i++) {
			double a = coefficients[0][i];
			double b = coefficients[1][i];
			double c = coefficients[2][i];
			double d = coefficients[3][i];
			//System.out.printf("(%8f, %8f, %8f, %8f)\n", a, b, c, d);
			System.out.printf("%12f + %12f(x - %2.0f) + %12f(x - %2.0f)^2 + %12f(x - %2.0f)^3\n", a, b, x[i], c, x[i], d, x[i]);
		}
		
		
		assertEquals(y[0], S.get(0).apply(x[0]), 1E-6);
		assertEquals(y[1], S.get(1).apply(x[1]), 1E-6);
		assertEquals(y[2], S.get(1).apply(x[2]), 1E-6);
		assertEquals(y[3], S.get(2).apply(x[3]), 1E-6);
		
	}
}
