package initial_value_ODE_methods;

import java.util.ArrayList;
import java.util.function.BiFunction;

public class DifferentialEquationSolver {
	
	public static ArrayList<Double> euler(double a, double b, int n, BiFunction<Double, Double, Double> f, double initValue) {
		double h = (b - a) / n;
		ArrayList<Double> w = new ArrayList<>();
		double w0 = initValue;
		w.add(w0);
		for (int i = 0; i < n; i++) {
			w0 = w0 + h * f.apply(a + i * h, w0);
			w.add(w0);
		}
		return w;
	}
	
	public static ArrayList<Double> modifiedEuler(double a, double b, int n, BiFunction<Double, Double, Double> f, double initValue) {
		double h = (b - a) / n;
		double w0 = initValue;
		ArrayList<Double> w = new ArrayList<>(n);
		w.add(w0);
		for (int i = 0; i <= n; i++) {
			double t = a + (i * h);
			w0 = w0 + (h / 2) * (f.apply(t, w0)+ f.apply(t + h, w0 + h * f.apply(t, w0)));
			w.add(w0);
		}
		return w;
	}
	
	public static ArrayList<Double> taylor(double a, double b, int n, BiFunction<Double, Double, Double> f, BiFunction<Double, Double, Double> fPrime, double initValue) {
		double h = (b - a) / n;
		double w0 = initValue;
		ArrayList<Double> w = new ArrayList<>();
		w.add(w0);
		/* 2nd order Taylor
		 * y(t[i+1]) = y(t[i]) + h*f(t[i], y(t[i])) + h^2/2 * f'(t[i], y(t[i])
		 * 
		 * w[i+1] = w[i] + hT^n(t[i],w[i]) for each i = 0 to N-1.
		 * T^n(t[i],w[i])= f(t[i],w[i]) + h/2 * f'(t[i],w[i]) + ... + (h^n-1 / n!) * f^(n-1)(t[i],w[i])
		 */
		for (int i = 0; i <= n; i++) {
			double t = a + (i * h);
			System.out.printf("t[i]= %.1f  %f\n", t, w0);
			//w[i+1] = w[i] + h * (f(t[i],w[i]) + h/2 * f'(t[i],w[i])
			w0 = w0 + h * (f.apply(t, w0) + (h/2 * fPrime.apply(t, w0)));
			w.add(w0);
		}
		return w;
	}
	
	public static ArrayList<Double> rungeKutta4(double a, double b, int n, BiFunction<Double, Double, Double> f, double initialValue) {
		ArrayList<Double> w = new ArrayList<>();
		double h = (b - a) / n;
		double w0 = initialValue;
		w.add(w0);
		for (int i = 0; i < n; i++) {
			double t = a + (i * h);
			double k1 = h * f.apply(t, w0);
			double k2 = h * f.apply(t + h / 2.0, w0 + (k1 / 2.0));
			double k3 = h * f.apply(t + h / 2.0, w0 + (k2 / 2.0));
			double k4 = h * f.apply(t + h, w0 + k3);
			w0 = w0 + (k1 + (2 * k2) + (2 * k3) + k4) / 6.0;
			w.add(w0);
		}
		return w;
	}
	
	public static ArrayList<Double> rungeKuttaFehlberg(double a, double b, BiFunction<Double, Double, Double> f, double initialValue, double tolerance, double hMax, double hMin) {
		ArrayList<Double> results = new ArrayList<>();
		
		double t = a;
		double w = initialValue;
		double h = hMax;
		boolean flag = true;
		
		results.add(w);
		System.out.printf("%12f  %12f\n", t, w);
		double k1, k2, k3, k4, k5, k6, r, e;
		int iterations = 0;
		while (flag) {
			iterations++;
			k1 = h * f.apply(t, w);
			k2 = h * f.apply(t + 0.25 * h, w + (k1 / 4.0));
			k3 = h * f.apply(t + ((3.0 / 8.0) * h), w + (3.0 / 32.0) * k1 + (9.0 / 32.0) * k2);
			k4 = h * f.apply(t + (12.0 / 13.0) * h, w + (1932.0 / 2197.0) * k1 - (7200.0 / 2197.0) * k2 + (7296.0 / 2197.0) * k3);
			k5 = h * f.apply(t + h, w + (439.0 / 216.0) * k1 - (8.0 * k2) + (3680.0 / 513.0) * k3 - (845.0 / 4104.0) * k4);
			k6 = h * f.apply(t + 0.5 * h, w - (8.0 / 27.0) * k1 + (2.0 * k2) - (3544.0 / 2565.0) * k3 + (1859.0 / 4104.0) * k4 - (11.0 / 40.0) * k5);
			
			r = (1.0 / h) * Math.abs((k1 / 360.0) - ((128.0 / 4275.0) * k3) - ((2197.0 / 75240.0) * k4) + ((1.0 / 50.0) * k5) + ((2.0 / 55.0) * k6));
			
			if (r <= tolerance) {
				t = t + h;
				w = w + (25.0 / 216.0) * k1 + (1408.0 / 2565.0) * k3 + (2197.0 / 4104.0) * k4 - (0.2) * k5;
				results.add(w);
				System.out.printf("%12f  %12f  %12f  %12.24f\n", t, w, h, r);
				
			}
			e = 0.84 * Math.pow(tolerance / r, 1.0 / 4.0);
			if (e <= 0.1) {
				h = 0.1 * h;
			}
			else if (e >= 4) {
				h = 4.0 * h;
			}
			else {
				h = e * h;
			}
			if (h > hMax) {
				h = hMax;
			}
			if (t >= b) {
				flag = false;
			}
			else if (t + h > b) {
				h = b - t;
			}
			else if (h < hMin) {
				flag = false;
			}
		}
		System.out.println(iterations);
		return results;
	}
}
































































