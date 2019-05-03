package differential_equation_tests;

import static org.junit.Assert.*;

import java.util.function.BiFunction;
import java.util.function.DoubleFunction;

import org.junit.Test;

import initial_value_ODE_methods.DifferentialEquationSolver;

public class EulerDiffTests {

	@Test
	public void testEulerDiff() {
		BiFunction<Double, Double, Double> f = (t, y) -> {
			return (y * y) / (1 + t);
		};
		double initialValue = - 1 / Math.log(2);
		double a = 1.0;
		double b = 2.0;
		int n = 10;
				
		DoubleFunction<Double> y = (t) -> { return - 1 / Math.log(t + 1); };
		
		DoubleFunction<Double> expected = (t) -> {
			return (y.apply(t) * y.apply(t)) / (1 + t);
		};
		
		displayDifferences(f, y, initialValue, a, b, n);
		DifferentialEquationSolver.euler(a, b, n, f, initialValue);
	}
	
	private static void displayDifferences(BiFunction<Double, Double, Double> f, DoubleFunction<Double> e, double init, double a, double b, int n) {
		double h = (b - a) / n;
		double w = init;
		for (int i = 0; i < n; i++) {
			double t = a + i * h;
			double actual = (double) e.apply(t);
			double error = actual - w;
			System.out.printf("%12f  -  %12f  =  %12f\n", w, actual, error);
			w = w + h * f.apply(t, w);
		}
	}
}
