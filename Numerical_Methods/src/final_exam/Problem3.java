package final_exam;

import java.util.ArrayList;
import java.util.function.BiFunction;
import java.util.function.DoubleFunction;

import initial_value_ODE_methods.DifferentialEquationSolver;

public class Problem3 {
	public static void main(String[] args) {
		part1();
		//part2();
	}
	
	private static void part1() {
		BiFunction<Double, Double, Double> yPrime = (t, y)-> {
			return ((y * y) + y) / t;
		};
		
		double h = 0.000001;
		double initialValue = -2.0;
		double a = 1.0;
		double b = 3.0;
		ArrayList<Double> results = DifferentialEquationSolver.rungeKutta4(a, b, (int) ((b - a) / h), yPrime, initialValue);
		for (int i = 0; i < results.size(); i++) {
			System.out.printf("%d. %.12f\n", i, results.get(i));
		}
		
		DoubleFunction<Double> exact = (t)-> {
			return (2 * t) / (1 - 2 * t);
		};
		
		System.out.printf("%.12f\n", exact.apply(b));
		System.out.printf("%.12f\n", Math.abs(exact.apply(b) - results.get(results.size() - 1)));
	}
	
	private static void part2() {
		BiFunction<Double, Double, Double> f = (t, y)-> {
			return ((y * y) + y) / t;
		};
		
		double initialValue = -2.0;
		double a = 1.0;
		double b = 3.0;
		double tolerance = Math.pow(10, -5);
		double hMax = 0.5;
		double hMin = 0.02;
		ArrayList<Double> results = DifferentialEquationSolver.rungeKuttaFehlberg(a, b, f, initialValue, tolerance, hMax, hMin);
		for (int i = 0; i < results.size(); i++) {
			System.out.printf("%d. %.12f\n", i, results.get(i));
		}
		
		DoubleFunction<Double> exact = (t)-> {
			return (2 * t) / (1 - 2 * t);
		};
		
		System.out.printf("Exact: %.12f\n", exact.apply(b));
		System.out.printf("Error: %.12f\n", Math.abs(exact.apply(b) - results.get(results.size() - 1)));
	}
}
