package final_exam;

import java.util.function.Function;

import numerical_integration.Integrate;

public class Problem1 {

	public static void main(String[] args) {
		Function<Double, Double> f = (x) -> {
			return (x * x) / (1 + Math.pow(x, 3));
		};
		
		double exact = Math.log(2) / 3;
		
		double a = 0.0;
		double b = 1.0;
		int n = 6;
		
		double[][] rombergResults = Integrate.romberg(a, b, f, n, 10E-7, exact);
		System.out.println(exact);
		System.out.println(rombergResults[2][n-1]);
		System.out.println(exact - rombergResults[2][n-1]);
		
		/*
		double simsponsResult = Integrate.compositeSimpson(a, b, n, f);
		System.out.println(exact);
		System.out.println(simsponsResult);
		System.out.println(exact - simsponsResult);
		*/
	}

}
