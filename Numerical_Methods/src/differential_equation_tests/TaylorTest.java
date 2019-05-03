package differential_equation_tests;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.function.BiFunction;
import java.util.function.DoubleFunction;

import org.junit.Test;

import initial_value_ODE_methods.DifferentialEquationSolver;

public class TaylorTest {

	@Test
	public void test() {
		BiFunction<Double, Double, Double> f = (t, y) -> {
			return (1 + t) / (1 + y);
		};
		BiFunction<Double, Double, Double> fPrime = (t, y) -> {
			return 1 / (y + 1);
		};
		double a = 1.0;
		double b = 2.0;
		double initValue = 2.0;
		int n = 10;
		ArrayList<Double> w = DifferentialEquationSolver.taylor(a, b, n, f, fPrime, initValue);
		
		DoubleFunction<Double> exact = (t)->{ 
			return Math.sqrt((t * t) + 2 * t + 6) - 1;
		};
		
		for (int i = 0; i <= n; i++) {
			double t = a + i * ((b - a) / 10);
			System.out.println(w.get(i) - exact.apply(t));
		}
	}
}
