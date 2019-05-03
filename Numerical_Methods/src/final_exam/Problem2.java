package final_exam;

import java.util.ArrayList;
import java.util.function.BiFunction;
import java.util.function.DoubleFunction;

import initial_value_ODE_methods.DifferentialEquationSolver;

public class Problem2 {

	public static void main(String[] args) {
		BiFunction<Double, Double, Double> p = (x, t) -> {
			return (.1) * (.02) * (1 - t);
		};
		ArrayList<Double> eulerValues = DifferentialEquationSolver.modifiedEuler(0, 50, 600, p, .01);
		for (int i = 0; i < eulerValues.size(); i++) {
			if (i % 12 == 0) {
				System.out.println(i + ". " + eulerValues.get(i));
			}
		}
		DoubleFunction<Double> exact = (t)-> {
			return 1 - (0.99 * Math.pow(Math.E, -.002 * t));
		};
		System.out.printf("Exact: %.12f\n", exact.apply(50));;
		System.out.printf("Error: %.12f\n", Math.abs(exact.apply(50) - eulerValues.get(eulerValues.size() - 2)));
	}
}