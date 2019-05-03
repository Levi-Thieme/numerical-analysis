package numerical_integration;

import java.util.function.Function;

/**
 * The Integrate class provides methods to approximate integrals of 
 * functions over specificed intervals. 
 * Various Newton-Cotes equations are utilized to approximate integrals.
 * 
 * @author Levi Thieme
 *
 */
public class Integrate {
	
	/**
	 * Applies the Trapezoidal rule to f on the interval [a, b] with n = 2;
	 * 
	 * @param a Left inclusive boundary.
	 * @param b Right inclusive boundary.
	 * @param f Double Function for which to approximate the integral over [a, b].
	 * @return The integral of f over [a, b].
	 */
	public static double trapezoid(double a, double b, Function<Double, Double> f) {
		double h = (b - a) / 2;
		return (h / 2) * (f.apply(a) + f.apply(b));
	}
	
	/**
	 * Calculates the error term of the trapezoidal rule
	 * 
	 * @param a Left inclusive boundary.
	 * @param b Right inclusive boundary.
	 * @param fPrime The second derivative of the function which has been integrated
	 * using the trapezoidal rule with n = 2 over [a, b].
	 * @param max The x-value such that fPrime(x) is maximum on the interval [a, b].
	 * @return The error term for the trapezoidal rule applied to a function f over the interval [a, b],
	 * where fPrime is the second derivative of f.
	 */
	public static double trapezoidError(double a, double b, Function<Double, Double> fPrime, double max) {
		double h = (b- a) / 2;
		return ((h * h * h) / 12) * fPrime.apply(max);
	}
	
	public static double compositeTrapezoid(double a, double b, int n, Function<Double, Double> f) {
		double h = (b - a) / n;
		double compositeSum = 0.0;
		for (int i = 1; i <= n - 1; i++) {
			compositeSum += f.apply(a + i * h);
		}
		return (h / 2) * (f.apply(a) + 2 * compositeSum + f.apply(b));
	}
	
	public static double compositeTrapezoidError(double a, double b, int n, Function<Double, Double> fPrime, double max) {
		double h = (b - a) / n;
		return (b - a) / 12 * (h * h) * fPrime.apply(max);
	}
	
	/**
	 * Integrates f over [a, b] using Simpson's rule with n = 2.
	 * 
	 * @param a Left inclusive boundary.
	 * @param b Right inclusive boundary.
	 * @param f Double function for which to approximate the integral over [a, b].
	 * @return The integral of f over [a, b].
	 */
	public static double simpson(double a, double b, Function<Double, Double> f) {
		double h = (b - a) / 2;
		return (h / 3) * (f.apply(a) + 4 * f.apply(a + h) + f.apply(b));
	}
	
	/**
	 * Computes the error for the Simpson's rule with n = 2.
	 * 
	 * @param a Left inclusive boundary.
	 * @param b Right inclusive boundary.
	 * @param fPrime The fourth derivative of the function which has been integrated
	 * using Simpson's rule with n = 2 over [a, b].
	 * @param max
	 * @return
	 */
	public static double simpsonError(double a, double b, Function<Double, Double> fPrime, double max) {
		double h = (b - a) / 2;
		return (Math.pow(h, 5) / 90) * fPrime.apply(max);
	}
	
	public static double simpsonThreeEighths(double a, double b, Function<Double, Double> f) {
		double h = (b - a) / 4;
		return ((3 * h) / 8) * (f.apply(a) + 3 * f.apply(a + h) + 3 * f.apply(a + 2 * h) + f.apply(b));
	}
	
	public static double simpsonsThreeEightsError(double a, double b, Function<Double, Double> fPrime, double max) {
		double h = (b - a) / 4;
		return ((3 * Math.pow(h, 5)) / 80) * fPrime.apply(max);
	}
	
	public static double closedNewtonCotesN4(double a, double b, Function<Double, Double> f) {
		double h = (b - a) / 4;
		return ((2 * h) / 45) * (7 * f.apply(a) + 32 * f.apply(a + h) + 12 * f.apply(a + 2 * h) + 32 * f.apply(a + 3 * h) + 7 * f.apply(b));
	}
	
	public static double closedNewtonCotesN4Error(double a, double b, Function<Double, Double> fPrime, double max) {
		double h = (b - a) / 4;
		return ((8 * Math.pow(h, 7)) / 945) * fPrime.apply(max);
	}
	
	public static double compositeSimpson(double a, double b, int n, Function<Double, Double> f) {
		double h = (b - a) / n;
		double compositeFastSum = 0.0;
		double compositeSlowSum = 0.0;
		
		for (int i = 1; i <= (n / 2) - 1; i++) {
			compositeFastSum += f.apply(a + h * (2 * i));
			compositeSlowSum += f.apply(a + h * ((2 * i) - 1));
		}
		
		compositeSlowSum += f.apply(a + h * ((2 * (n / 2) - 1)));
		
		return (h / 3) * (f.apply(a) + 2 * compositeFastSum + 4 * compositeSlowSum + f.apply(b));
	}
	
	public static double compositeSimpsonError(double a, double b, int n, Function<Double, Double> fPrime, double max) {
		double h = (b - a) /  n;
		return (b - a) / 180 * Math.pow(h, 4) * fPrime.apply(max);
	}
	
	public static double[][] romberg(double a, double b, Function<Double, Double> f, int n, double tolerance, double exact) {
		double[][] r = new double[3][n];
		double h = b - a;
		r[1][1] = (h / 2.0) * (f.apply(a) + f.apply(b));
		System.out.println(r[1][1]);
		for (int i = 2; i < n; i++) {
			double sum = 0.0;
			for (double k = 1; k <= Math.pow(2, i - 2); k++) {
				sum += f.apply(a + ((k - 0.5) * h));
			}
			
			r[2][1] = 0.5 * (r[1][1] + (h * sum));
			
			for (int j = 2; j <= i; j++) {
				r[2][j] = r[2][j-1] + ((r[2][j-1] - r[1][j-1]) / (Math.pow(4.0, j - 1) - 1));
			}
			for (int j = 1; j <= i; j++) {
				System.out.printf("%16.10f", r[2][j]);
			}
			System.out.println("\n");
			h = h / 2.0;
			for (int j = 1; j < i; j++) {
				r[1][j] = r[2][j];
			}
			if (Math.abs(r[1][i]- exact) < tolerance) {
				break;
			}
		}
		return r;
	}
}





































