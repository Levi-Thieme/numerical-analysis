package tests;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.function.BiFunction;
import java.util.function.DoubleFunction;

import org.junit.Test;

import initial_value_ODE_methods.DifferentialEquationSolver;

public class RFKTest {

	@Test
	public void test() {
		BiFunction<Double, Double, Double> yPrime = (t, y)-> y - (t * t) + 1;
		DoubleFunction<Double> exact = (t)-> Math.pow(t + 1, 2) - (0.5 * Math.pow(Math.E, t));
		double tolerance = Math.pow(10, -5);
		double hMax = 0.25;
		double hMin = 0.01;
		double initialValue = 0.5;
		double a = 0.0;
		double b = 2.0;
		
		ArrayList<Double> results = DifferentialEquationSolver.rungeKuttaFehlberg(a, b, yPrime, initialValue, tolerance, hMax, hMin);
		double[] expected = { 0.5, 0.9204886, 1.3964910, 1.9537488, 2.5864260, 3.2604605,
							  3.9520955, 4.6308268, 5.2574861, 5.3054896 };
		for (int i = 0; i < results.size(); i++) {
			assertEquals(expected[i], results.get(i), 1E-7);
		}
	}

}
