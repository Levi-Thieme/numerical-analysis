package numerical_functions;

import java.util.ArrayList;
import java.util.function.BiFunction;
import java.util.function.Function;

public class NumericalFunction {
	/*pseudocode for calculating a^pow mod b
	 * 
	 * if (power == 1) then return a mod b;
	 * else if (power odd) then return a * this(a, pow - 1);
	 * else then return this(a, pow / 2) * this(a, pow / 2);
	 */
	
	/**
	 * Calculates the modulus of an exponential constant base^pow mod b.
	 * 
	 * @param base The base of the exponential constant term.
	 * @param pow The power of the exponential constant term.
	 * @param mod The modulus value.
	 * @return (base^pow) % mod
	 * O(log(n))
	 */
	public static double exponentialModulus(double base, int pow, int mod) {
		if (pow == 1) { return base % mod; }
		else if (pow % 2 == 1) { return ((base % mod ) * (exponentialModulus(base, pow - 1, mod)) % mod); }
		int powerHalf = pow / 2;
		return ((exponentialModulus(base, powerHalf, mod) % mod) * (exponentialModulus(base, powerHalf, mod) % mod)) % mod;
	}
	
	/**
	 * Calculates the modulus of an exponential constant x(^y^z) % p
	 * 
	 * @param x The base.
	 * @param y The exponent of base.
	 * @param z The exponent of y.
	 * @param p The modulus value.
	 * @return x^(y^z) % p
	 * O(z * nlog(n))
	 */
	public static double exponentialExponentialModulus(double x, int y, int z, int p) {
		int powerDelta = 1;
		// Iteratively apply (x^a % p) ^ y^(z-a) % p where a is powerDelta
		while (z > 1) {
			z = z - powerDelta;
			x = exponentialModulus(x, (int)Math.pow(y, powerDelta), p);
			//this results in the calculation of x^a % p
		}
		return exponentialModulus(x, y, p);
	}
	
	/**
	 * Computes the greatest common divisor of a and b using Euclid's algorithm.
	 * 
	 * @param a The first dividend.
	 * @param b The second dividend.
	 * @return The greatest common divisor of a and b.
	 */
	public static long gcd(long a, long b) {
		if (a == 0) { return b; }
		if (b == 0) { return a; }
		//Rewrite a = b * q + c
		long r = a % b;
		return gcd(b, r);
	}
	
	/**
	 * The extended version of Euclid's gcd algorithm.
	 * ax + by = gcd(a, b)
	 * 
	 * @param a The first dividend.
	 * @param b The second dividend.
	 * @param x 
	 * @param y
	 * @return The gcd(a, b) and integer coefficients a and b
	 * such that ax + by = gcd(a, b)
	 */
	public static long gcdExtended(long a, long b, long[] x) {
		if (a == 0) {
			x[0] = 0;
			x[1] = 1;
			return b;
		}
		long[] t = {1, 1};
		long gcd = gcdExtended(b % a, a, t);
		x[0] = t[1] - (b / a) * t[0];
		x[1] = t[0];
		return gcd;
	}
	
	/**
	 * Calculates the modular multiplicative inverse
	 * 
	 * @param a coefficient of a gcd factor
	 * @param m modulo value
	 * @return x such that a * x congruent 1 mod m
	 */
	public static long modInverse(long a, long m) {
		long[] t = {1, 1};
		long g = gcdExtended(a, m, t);
		if (g != 1) {
			System.out.println("Inverse does not exist.");
			return -1;
		}
		long inverse = (t[0] % m + m) % m;
		return inverse;
	}
	
	/**
	 * Returns a function(a degree 1 polynomial) for the interpolation between two points.
	 * 
	 * @param xVals An Array containing x0 and x1
	 * @param yVals An array containing y0 and y1
	 * @return A function that interpolates the xVals and yVals respectively.
	 */
	public static Function<Double, Double> firstDegreeInterp(double[] xVals, double[] yVals) {
		Function<Double, Double> f = (xInput) -> {
			double x0 = xVals[0];
			double x1 = xVals[1];
			double y0 = yVals[0];
			double y1 = yVals[1];
			double slope = (y1 - y0) / (x1 - x0);
			return slope * (xInput - x0) + y0;
		};
		return f;
	}
	
	/** TODO make this general, so as to return a function to be evaluated at x
	 * Computes the Lagrange interpolation between two points evaluated at x
	 * 
	 * P(x) = L0(x)f(x0) + L1(x)f(x1)
	 * 
	 * @param x Input for P
	 * @param xVals x0, x1
	 * @param yVals y0, y1
	 * @return P(x)
	 */
	public static double lagrangeInterpFirstDegree(double x, double[] xVals, double[] yVals) {
		/*
		BiFunction<Double, Double[], Double> L0 = (inp, xInputs) -> {
			return (x - xVals[1]) / (xVals[0] - xVals[1]);
		};
		
		BiFunction<Double, Double[], Double> L1 = (inp, xInputs) -> {
			return (x - xVals[0]) / (xVals[1] - xVals[0]);
		};
		*/
		double L0 = (x - xVals[1]) / (xVals[0] - xVals[1]);
		double L1 = (x - xVals[0]) / (xVals[1] - xVals[0]);
		Function<Double, Double> f = firstDegreeInterp(xVals, yVals);
		//L0(x)f(x0) + L1(x)f(x1)
		return (L0 * f.apply(xVals[0])) + (L1 * f.apply(xVals[1]));
	}
	
	
	/**
	 * The lagrangeInterp function returns a function that interpolates n (x, y) pairs.
	 * This implementation is the "naive" approach, and probably isn't suitable for large degrees. 
	 * Barycentric interpolation provides a more time efficient algorithm for interpolation.
	 * 
	 * 
	 * sum(0...n) of f(x[k]) * L[k](x)
	 * 
	 * @param x A double array containing n x-values.
	 * @param y A double array containing n y-values.
	 * @param coefficientPolynomials An ArrayList of Bifunctions that compute the polynomial coefficients
	 * @return A function that computes sum(0...n) of f(x[k]) * L[k](x)
	 */
	public static Function<Double, Double> lagrangeInterp(Double[] xVals, double[] yVals, ArrayList<BiFunction<Double[], Double[], Double>> coefficientPolynomials) {
		Function<Double, Double> f = (x) -> {
			double total = 0.0;
			for (int i = 0; i < xVals.length; i++) {
				Double[] xAndK = {x, (double) i};
				total += yVals[i] * coefficientPolynomials.get(i).apply(xVals, xAndK);
			}
			return total;
		};
		return f;
	}
	
	/**
	 * This function returns an ArrayList containing the coefficient polynomials for Lagrange polynomial interpolation.
	 * The coefficient polynomials are the L functions in the definition of a Lagrange polynomial.
	 * Lagrange polynomial: sum of n polynomials where each polynomial is f(x[k]) * L(x)
	 * 
	 * @param xVals A double array containing the x values.
	 * @param n The degree of the coefficient polynomials.
	 * @return An ArrayList containing the coefficient polynomials.
	 */
	public static ArrayList<BiFunction<Double[], Double[], Double>> lagrangeCoefficientPolynomials(int n) {
		ArrayList<BiFunction<Double[], Double[], Double>> coefficientPolynomials = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			/* f is the kth Lagrange coefficient polynomial and is
			 * defined as the following: for every j != k compute (x - x[j]) / (x[k] - x[j])
			 */
			BiFunction<Double[], Double[], Double> f = (xVals, xAndK) -> {
				double total = 1.0;
				double x = xAndK[0];
				double k = xAndK[1];
				for (int j = 0; j < n; j++) {
					if (j != k) {
						total *= (x - xVals[j]) / (xVals[(int)k] - xVals[j]);
					}
				}
				return total;
			};
			coefficientPolynomials.add(f);
		}
		return coefficientPolynomials;
	}
	
	public static double lagrangeInterpError(int a, int b, Function<Double, Double> nthDerivative, Function<Double, Double> polynomialDerivative) throws Exception {
		//find max value of nthDerivative, test endpoints and where the derivative is zero, get the max of these 3 values
		
		//compute the absolute maximum value of polynomialDerivative
		throw new Exception("Not implemented yet.");
	}
	
	/**
	 * Calculates the output of f given a list of x-values and returns the outputs.
	 * 
	 * @param x A double array of x-values
	 * @param f The function to be evaluated for all x-values
	 * @return An array containing the values of f evaluated at every x-value
	 */
	public static double[] evaluateFunction(double[] x, Function<Double, Double> f) {
		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			y[i] = f.apply(x[i]);
		}
		return y;
	}
	
	/**
	 * Returns the divided differences of x[i]...x[j]
	 * 
	 * @param x A double array of x-values
	 * @param f The function to be evaluated
	 * @return j - i divided differences
	 */
	public static double[][] dividedDiff(Double[] x,  Double[] y) {
		int n = x.length;
		double[][] differences = new double[n][n];
		//Fill differences matrix first column with f(X0,f(X1),...,f(Xn)
		for (int i = 0; i < y.length; i++) {
			differences[i][0] = y[i];
		}
		/*Compute the divided differences
		 * 
		 * F[i][j] = (F[i][j-1] - F[i-1][j-1]) / (x[i] - x[i-j])
		 */
		for (int i = 1; i < n; i++) {
			for (int j = 1; j <= i; j++) {
				differences[i][j] = (differences[i][j-1] - differences[i-1][j-1]) / (x[i] - x[i-j]);
			}
		}
		return differences;
	}
	
	/**
	 * Returns the divided differences of x[i]...x[j]
	 * Precondition: differences must contain f(X0),f(X1),..,f(Xn) in first column.
	 * @param x A double array of x-values
	 * @param f The function to be evaluated
	 * @return j - i divided differences
	 */
	public static double[][] dividedDiff(double[] x, double[][] differences) {
		int n = x.length;
				
		for (int i = 1; i < n; i++) {
			for (int j = 1; j <= i; j++) {
				differences[i][j] = (differences[i][j-1] - differences[i-1][j-1]) / (x[i] - x[i-j]);
			}
		}
		return differences;
	}
	
	/**
	 * Creates a divided differences matrix with the next divided differences.
	 * 
	 * @param x The x-values.
	 * @param differences Matrix containing the divided differences.
	 * @return A double matrix containing the differences and the next divided difference.
	 */
	public static double[][] extendDividedDifferences(Double[] x, Double[] newY, double[][] differences) {
		
		/*
		 * Construct a new row of divided differences. F[n+1][0]F[n+1][1],...,F[n+1][n+1]
		 */
		int maxRow = differences.length + 1;
		int maxCol = differences[differences.length-1].length + 1;
		double[][] diffs = new double[maxRow][maxCol];
		//copy differences to larger matrix
		for (int i = 0; i < maxRow - 1; i++) {
			for (int j = 0; j < maxCol - 1; j++) {
				diffs[i][j] = differences[i][j];
			}
		}
		//compute the new row of divided differences F[n+1][0]...F[n+1][n+1]
		int row = maxRow - 1;
		diffs[row][0] = x[x.length-1];
		for (int i = 1; i < maxCol; i++) {
			diffs[row][i] = ((i-1 == 0 ? newY[row] : diffs[row][i-1]) - diffs[row-1][i-1]) / (x[row] - x[row-i]);
		}
		return diffs;
	}
	
	/**
	 * 
	 * @param x x-values for a function f
	 * @param y y-values for a function f
	 * @param i High index
	 * @param j Low index
	 * @return The divided difference of x and y for i and j.
	 */
	public static double dividedDiff(double[] x, double[] y, int i, int j) {
		return (y[j] - y[i]) / (x[j] - x[i]);
	}
	
	
	public static double lagrangianPolynomialDenominator(Double x, Double[] xVals, int n) {
		double total = 1.0;
		for (int i = 0; i < n; i++) {
			System.out.printf("(%12f - %12f)", x, xVals[i]);
			total *= (x - xVals[i]);
		}
		System.out.println("\n");
		return total;
	}
	
	/**
	 * Interpolates an array of (x,y) pairs using Newton's divided differences.
	 * 
	 * @param x Double array of x-values.
	 * @param y Double array of y-values.
	 * @return A Function<Double, Double> that interpolates the (x,y) values.
	 */
	public static Function<Double, Double> newtonInterp(Double[] x, Double[] y) {
		//Generate the divided differences
		final double[][] differences = NumericalFunction.dividedDiff(x, y);
		//Pn(X) = F00 + sum (Fii product(x-xj))
		final int n = x.length;
				
		Function<Double, Double> f = (xVal) -> {
			Double total = differences[0][0];
			for (int i = 1; i < n; i++) {
				double coefficientPoly = NumericalFunction.lagrangianPolynomialDenominator(xVal, x, i);
				total += differences[i][i] * coefficientPoly;
				System.out.printf("%d.  %f * %f\n", i, differences[i][i], coefficientPoly);
			}
			return total;
		};
		return f;
	}
	
	/**
	 * Interpolates an array of (x,y) pairs using Newton's divided differences.
	 * 
	 * @param x Double array of x-values.
	 * @param y Double array of y-values.
	 * @param coefficientPolynomials An ArrayList of the coefficient polynomials: (x-x0)(x-x1)...(x-xk)
	 * @param dividedDifferences A double matrix whose diagonal contains the divided differences of the Double array x.
	 * @return A Function<Double, Double> that interpolates the (x,y) values.
	 */
	public static Function<Double, Double> newtonInterp(Double[] x, Double[] y, ArrayList<BiFunction<Double[], Double[], Double>> coefficientPolynomials, double[][] dividedDifferences) {
		Function<Double, Double> f = (xVal) -> {
			Double total = 0.0;
			for (int i = 0; i < dividedDifferences[0].length; i++) {
				double diff = dividedDifferences[i][i];
				double coefficientPoly = NumericalFunction.lagrangianPolynomialDenominator(xVal, x, i);
				total += diff * coefficientPoly;
			}
			return total;
		};
		return f;
	}
	
	/**
	 * Attempts to find a fixed point for the function f within the interval [a,b].
	 * 
	 * @param f A continuous function over the the interval [a,b]
	 * @param a Left interval endpoint
	 * @param b Right interval endpoint
	 * @param n The maximum number of iterations
	 * @param x0 The starting value for the fixed point iteration
	 * @return An arraylist containing the fixed points calculated at each iteration.
	 */
	public static ArrayList<Double> fixedPoint(Function<Double, Double> f, double a, double b, int n, double x0, double tolerance) {
		ArrayList<Double> points = new ArrayList<>();
		int i = 0;
		double p = x0;
		points.add(p);
		double p0 = Double.MAX_VALUE;
		while (i < n) {
			p0 = f.apply(p);
			points.add(p0);
			i += 1;
			if (Math.abs(p0 - p) <  tolerance) {
				break;
			}
			p = p0;
		}
		return points;
	}
	
	/**
	 * Returns a matrix containing the coefficients for the natural cubic spline interpolation of the (x,y) values.
	 * 
	 * @param x A Double array of x-values.
	 * @param a A Double array of y-values.
	 * @return A matrix containing the coefficients a,b,c,d for indices 0,1,...,n-1
	 */
	public static double[][] naturalCubicSplineCoefficients(Double[] x, double[] a) {
		int n = x.length;
		//Compute h0...hj
		double[] h = new double[n];
		for (int i = 0; i < n - 1; i++) {
			h[i] = x[i+1] - x[i];
		}
		
		//Compute a0...aj
		double[] alpha = new double[n];
		alpha[n-1] = 0.0;
		for (int i = 1; i < n - 1; i++) {
			alpha[i] = ((3 / h[i]) * (a[i+1] - a[i])) - ((3 / h[i-1]) * (a[i] - a[i-1]));
		}
		
		
		//Solve the tridiagonal linear system for constants
		double[] l = new double[n];
		double[] u = new double[n];
		double[] z = new double[n];
		l[0] = 1.0;
		u[0] = 0.0;
		z[0] = 0.0;
		
		for (int i = 1; i < n - 1; i++) {
			l[i] = (2 * (x[i+1] - x[i-1])) - (h[i-1] * u[i-1]);
			u[i] = h[i] / l[i];
			z[i] = (alpha[i] - (h[i-1] * z[i-1])) / l[i];
		}
		
		double[] b = new double[n];
		double[] c = new double[n];
		double[] d = new double[n];
		l[n-1] = 1.0;
		z[n-1] = 0.0;
		c[n-1] = 0.0;
		
		for (int j = n-2; j >= 0; j--) {
			c[j] = z[j] - (u[j] * c[j+1]);
			b[j] = ((a[j+1] - a[j]) / h[j]) - (h[j] * (c[j+1] + (2 * c[j]))) / 3;
			d[j] = (c[j+1] - c[j]) / (3 * h[j]);
		}
		
		double[][] coefficients = new double[4][n];
		coefficients[0] = a;
		coefficients[1] = b;
		coefficients[2] = c;
		coefficients[3] = d;
		return coefficients;
	}
	
	/**
	 * Computes the coefficients, a,b,c,d for 0,1,...,n, for a clamped cubic spline interpolation polynomial
	 * 
	 * @param x X-values
	 * @param a Y-values
	 * @param fp0 = f'(x[0])
	 * @param fpn = f'(x[n])
	 * @return A matrix containing the constants a,b,c,d
	 */
	public static double[][] clampedCubicSplineCoefficients(Double[] x, double[] a, double fp0, double fpn) {
		int n = x.length;
		
		//Compute h-values
		double[] h = new double[n];
		for (int i = 0; i < n - 1; i++) {
			h[i] = x[i+1] - x[i];
		}
		
		//Compute a-values
		double[] alpha = new double[n];
		alpha[0] = (3 * (a[1] - a[0]) / h[0]) - (3 * fp0);
		alpha[n-1] = (3 * fpn) - (3 * (a[n-1] - a[n-2]) / h[n-2]);
		
		for (int i = 1; i < n - 1; i++) {
			alpha[i] = (3 / h[i]) * (a[i+1] - a[i]) - (3 / h[i-1]) * (a[i] - a[i-1]);
		}
		
		double[] l = new double[n];
		double[] u = new double[n];
		double[] z = new double[n];
		double[] b = new double[n];
		double[] c = new double[n];
		double[] d = new double[n];
		
		l[0] = 2 * h[0];
		u[0] = 0.5;
		z[0] = alpha[0] / l[0];
		
		for (int i = 1; i < n-1; i++) {
			l[i] = 2 * (x[i+1] - x[i-1]) - (h[i-1] * u[i-1]);
			u[i] = h[i] / l[i];
			z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
		}
		
		l[n-1] = h[n-2] * (2 - u[n-2]);
		z[n-1] = (alpha[n-1] - h[n-2] * z[n-2]) / l[n-1];
		c[n-1] = z[n-1];
		
		for (int j = n-2; j > 0; j--) {
			c[j] = z[j] - u[j] * c[j+1];
			b[j] = ((a[j+1] - a[j]) / h[j]) - (h[j] * (c[j+1] + 2 * c[j]) / 3);
			d[j] = (c[j+1] - c[j]) / (3 * h[j]);
		}
		
		double[][] coefficients = new double[4][];
		coefficients[0] = a;
		coefficients[1] = b;
		coefficients[2] = c;
		coefficients[3] = d;
		
		return coefficients;
	}
	
	
	/**
	 * Constructs an ArrayList of functions composed of cubic splines.
	 * 
	 * @param xVals Array of x-values
	 * @param coefficients a,b,c,d coefficients
	 * @param j Number of cubic splines
	 * @return An ArrayList of cubic functions that interpolate xValues
	 */
	public static ArrayList<Function<Double, Double>> cubicSplineInterp(Double[] xVals, double[][] coefficients, int j) {
		double[] a = coefficients[0];
		double[] b = coefficients[1];
		double[] c = coefficients[2];
		double[] d = coefficients[3];
		
		ArrayList<Function<Double, Double>> s = new ArrayList<>(j);
		for (int i = 0; i < j; i++) {
			final int k = i;
			Function<Double, Double> f = (x) -> {
				double h = x - xVals[k];
				return a[k] + b[k] * h + c[k] * (h * h) + d[k] * (h * h * h);
			};
			s.add(f);
		}
		return s;
	}
	
	public static double threePointEndpoint(Function<Double, Double> f, double x0, double h) {
		return (1 / (2 *h)) * (-3 * f.apply(x0) + 4 * f.apply(x0 + h) - f.apply(x0 + 2 * h));
	}
	
	public static double threePointMidpoint(Function<Double, Double> f, double x0, double h) {
		return (1 / (2 * h)) * (f.apply(x0 + h) - f.apply(x0 - h));
	}
	
	public static double fivePointMidpoint(Function<Double, Double> f, double x0, double h) {
		return (1 / (12 * h)) * (f.apply(x0 - (2 * h)) - (8 * f.apply(x0 - h)) + (8 * f.apply(x0 + h)) - f.apply(x0 + (2 * h)));
	}
	
	public static double fivePointEndpoint(Function<Double, Double> f, double x0, double h) {
		return (1 / (12 * h)) * (-25 * f.apply(x0) + 48 * f.apply(x0 + h) - 36 * f.apply(x0 + (2 * h)) + 
				16 * f.apply(x0 + (3 * h)) - 3 * f.apply(x0 + (4 * h)));
	}
}

















































