package exam1problems;

import java.util.ArrayList;
import java.util.function.Function;

import numerical_functions.NumericalFunction;

public class Problem1 {
	//Fixed point iteration
	public static void main(String[] args) {
		
		double x0 = 2.0;
		Function<Double, Double> f = (x)-> { return ((24 * x) + ( 25 / (x * x))) / 25; };
		double a = 1;
		double b = 3;
		double tolerance = 1E-12;
		int n = 500;
		/*
		ArrayList<Double> points = NumericalFunction.fixedPoint(f, a, b, n, x0, tolerance);
		System.out.println("Iterations " + points.size());
		printPoints(points);
		*/
	 
		Function<Double, Double> g = (x)-> { return ((2 * x * x * x) + 25) / ( 3 * x * x); };
		ArrayList<Double> points = NumericalFunction.fixedPoint(f, a, b, n, x0, tolerance);
		printPoints(points);
	}	
	
	private static void printPoints(ArrayList<Double> points) {
		int stepCount = 1;
		for (Double point : points) {
			System.out.printf("%3d %14f\n", stepCount,  point);
			stepCount++;
		}
	}
}
