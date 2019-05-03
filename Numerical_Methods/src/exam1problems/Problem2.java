package exam1problems;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.function.Function;

import root_finding_algorithms.RootCalculator;

public class Problem2 {

	public static void main(String[] args) {
		Function<Number, Number> f = (xValue)->{
			double x = xValue.doubleValue();
			return (x*x*x) + (2*x*x) + (10*x) - 20;
		};
		
		Function<Number, Number> fPrime = (xValue)-> {
			double x = xValue.doubleValue();
			return (3 * x * x) + (4 * x) + 10;
		};
		
		HashMap<String, Object> results = RootCalculator.newton(f, fPrime, 1, 500, 1E-16);
		ArrayList<Double> intercepts = (ArrayList<Double>) results.get("xIntercepts");
		if (intercepts != null) {
			for (Double x : intercepts) {
				System.out.println(x);
			}
		}
	}

}
