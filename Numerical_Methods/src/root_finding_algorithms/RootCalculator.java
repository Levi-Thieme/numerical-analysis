package root_finding_algorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.function.Function;

public class RootCalculator {
	
	/**
	 * Finds the root of a funciton using the Bisection algorithm.
	 * 
	 * @param f - the function for which to find a root
	 * @param a - inclusive lower interval bound
	 * @param b - inclusive upper interval bound
	 * @param maxSteps - the maximum number of iterations
	 * @param accuracy - an upper bound for the accuracy of the interval containing the root
	 * @return A HashMap containing the results of the algorithm
	 */
	public static HashMap<String, Object> bisection(Function<Number, Number> f, double a, double b, int maxSteps, double accuracy) {
		//print solution, step count, table of midpoints
		HashMap<String, Object> returnVals = new HashMap<>();
		ArrayList<Double> midpoints = new ArrayList<>();
		double fA = f.apply(a).doubleValue();
		double fB = f.apply(b).doubleValue();
		//Check that f(a) and f(b) have opposite signs
		if (fA * fB > 0) {
			return null;
		}
		
		int n = 1;
		double midpoint = 0.0;
		while (n <= maxSteps) {
			double half = (b - a) / 2.0;
			midpoint = a + half;
			midpoints.add(midpoint);
			double fM = f.apply(midpoint).doubleValue();
			//If we have found an exact root or our root interval is within accuracy
			if (fM == 0 || half < accuracy) {
				returnVals.put("accuracy", (Double)(accuracy - half));
				break;
			}
			n += 1;
			//Choose upper or lower half for the n+1 iteration
			if (fA * fM > 0) {
				a = midpoint;
			}
			else {
				b = midpoint;
			}
		}
		
		if (n > maxSteps) {
			returnVals.put("message", "Failed to find a root within " + maxSteps + " iterations.\n");
			returnVals.put("rootFound", false);
		}
		else {
			returnVals.put("message", "Found the root, " + midpoint + ", at " + n + " iterations with an accuracy of " + returnVals.get("accuracy"));
			returnVals.put("rootFound", true);
		}
		returnVals.put("steps", n);
		returnVals.put("root", midpoint);
		returnVals.put("midpoints", midpoints);
		return returnVals;
	}
	
	/**
	 * Prints the results of the RootBisection.findRoot function to the console.
	 * 
	 * @param results HashMap containing the results.
	 * @param fOfRoot The function evaluated at the approximated root.
	 */
	public static void displayRootBisectionResults(HashMap<String, Object> results, double fOfRoot) {
		if (results == null) {
			System.out.println("Invalid function interval. The function must have opposite signs at each boundary of the interval.");
			System.exit(0);
		}
		System.out.println(results.get("message"));
		if ((boolean) results.get("rootFound")) {
			System.out.println("f(root) = " + fOfRoot);
		}
		if (results.get("midpoints") instanceof ArrayList) {
			ArrayList<?> midpoints = (ArrayList<?>)results.get("midpoints");
			StringBuilder stringBuilder = new StringBuilder();
			for(int i = 0; i < midpoints.size(); i++) {
				stringBuilder.append(i + 1 + ".  " + midpoints.get(i) + "\n");
			}
			System.out.println(stringBuilder.toString());
		}
	}


	/**
	 * Finds the root of a function using Newton's algorithm.
	 * 
	 * @param f - the function for which to find a root
	 * @param fPrime - the first derivative of the function f
	 * @param start - the initial value to start the iteration at
	 * @param maxSteps - the maximum number of iterations
	 * @param accuracy - an upper bound for the accuracy of the interval containing the root
	 * @return A HashMap containing the results of the algorithm
	 */
	public static HashMap<String, Object> newton(Function<Number, Number> f, Function<Number, Number> fPrime, double start, int maxSteps, double accuracy) {
		// p1 = f(p0) - f(p0)/f'(p0)
		// p1 - po < accuracy then done
		HashMap<String, Object> results = new HashMap<>();
		ArrayList<Double> tangentXIntercepts = new ArrayList<>();
		int iterations = 0;
		double current= start;
		while (iterations < maxSteps) {
			double next = current - (f.apply(current).doubleValue() / fPrime.apply(current).doubleValue());
			if (Math.abs(current - next) < accuracy) {
				results.put("root", next);
				tangentXIntercepts.add(next);
				break;
			}
			current = next;
			tangentXIntercepts.add(current);
			iterations += 1;
		}
		results.put("xIntercepts", tangentXIntercepts);
		return results;
	}
	
	/**
	 * Finds the root of a function using the secant method.
	 * 
	 * @param f The function for which to find a root
	 * @param p0 The first root approximation
	 * @param p1 The second root approximation
	 * @param maxSteps The maximum number of iterations.
	 * @param accuracy The required accuracy for the root's approximation.
	 * @return A list of roots.
	 */
	public static ArrayList<Double> secant(Function<Number, Number> f, double p0, double p1, int maxSteps, double accuracy) {
		ArrayList<Double> roots = new ArrayList<>();
		double q0 = f.apply(p0).doubleValue();
		double q1 = f.apply(p1).doubleValue();
		int i = 2;
		while (i <= maxSteps) {
			double p = p1 - q1 * (p1 - p0) / (q1 -q0);
			roots.add(p);
			if (Math.abs(p - p1) < accuracy) {
				break;
			}
			i += 1;
			p0 = p1;
			q0 = q1;
			p1 = p;
			q1 = f.apply(p).doubleValue();
		}
		return roots;
	}
}























































