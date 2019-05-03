package differential_equation_tests;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.function.BiFunction;
import java.util.function.DoubleFunction;

import org.junit.Test;

import initial_value_ODE_methods.DifferentialEquationSolver;

public class RungeKutta4Test {

	@Test
	public void test() {
		double a = 1.0;
		double b = 2.0;
		int n = 2;
		double initValue = 2.0;
		
		BiFunction<Double, Double, Double> f = (t, y) -> {
			return (1 + t) / (1 + y);
		};
		
		ArrayList<Double> w = DifferentialEquationSolver.rungeKutta4(a, b, n, f, initValue);
		
		ArrayList<Double> expected = new ArrayList<>();
		DoubleFunction<Double> exact = (t)->{ 
			return Math.sqrt((t * t) + 2 * t + 6) - 1;
		};
		expected.add(initValue);
		for (int i = 1; i <= n; i++) {
			double t = a + i * ((b - a) / n);
			expected.add(exact.apply(t));
		}
		printError(n, w, expected);
	}
	
	private static void printError(int n, ArrayList<Double> actual, ArrayList<Double> expected) {
		for (int i = 0; i <= n; i++) {
			System.out.printf("%f - %f = %f\n", expected.get(i), actual.get(i), expected.get(i) - actual.get(i));
		}
	}
}
