package exam1problems;

import java.util.function.Function;

import numerical_functions.NumericalFunction;

public class Problem5 {

	public static void main(String[] args) {
		partA();
	}
	
	private static void partA() {
		Function<Double, Double> f = (x) -> {
			return Math.tan(2 * x);
		};
		
		//2sec^2(2x)
		Function<Double, Double> fPrime = (x) -> {
			return 2 * Math.pow(( 1 /Math.cos(2 * x)), 2);
		};
		
		double h = 0.05;
		double x0 = 0.0;
		
		/*
		//5 point midpoint to x= 1.15, 1.20
		double fivePointMid0 = NumericalFunction.fivePointMidpoint(f, 1.15, h);
		double fivePointMid1 = NumericalFunction.fivePointMidpoint(f, 1.20, h);
		System.out.printf("5p mid x= 1.15, y =%f\n5p mid x= 1.20, y =%f\n", fivePointMid0, fivePointMid1);
		
		//3 point midpoint to x=1.10,1.15,1.20,1.25
		double threePointMid = NumericalFunction.threePointMidpoint(f, 1.10, h);
		System.out.println(threePointMid);
		threePointMid = NumericalFunction.threePointMidpoint(f, 1.15, h);
		System.out.println(threePointMid);
		threePointMid = NumericalFunction.threePointMidpoint(f, 1.20, h);
		System.out.println(threePointMid);
		threePointMid = NumericalFunction.threePointMidpoint(f, 1.25, h);
		System.out.println(threePointMid);
		
		//3 point endpoint to x= 1.05, 1.30
		double threePointEnd = NumericalFunction.threePointEndpoint(f, 1.05, h);
		System.out.println(threePointEnd);
		threePointEnd = NumericalFunction.threePointEndpoint(f, 1.30, -h);
		System.out.println(threePointEnd);
		*/
		//5 point endpoint to x = 1.05, 1.10, 1.25, 1.30
		x0 = 1.05;
		double fivePointEnd = NumericalFunction.fivePointEndpoint(f, x0, h);
		p(fivePointEnd);
		p(error(fPrime.apply(x0), fivePointEnd));
		
		x0 = 1.10;
		fivePointEnd = NumericalFunction.fivePointEndpoint(f, x0, h);
		p(fivePointEnd);
		p(error(fPrime.apply(x0), fivePointEnd));
		
		x0 = 1.25;
		fivePointEnd = NumericalFunction.fivePointEndpoint(f, x0, -h);
		p(fivePointEnd);
		p(error(fPrime.apply(x0), fivePointEnd));
		
		x0 = 1.30;
		fivePointEnd = NumericalFunction.fivePointEndpoint(f, x0, -h);
		p(fivePointEnd);
		p(error(fPrime.apply(x0), fivePointEnd));
	}
	
	private static void partB() {
		Function<Double, Double> f = (x) -> {
			return (1 / Math.pow(Math.E, x)) - 1 + x;
		};
		
		Function<Double, Double> fPrime = (x) -> {
			return (-1 / Math.pow(Math.E, x)) + 1;
		};
		
		double h = 0.2;
		double x0 = 0.0;
		/*
		//5 point midpoint x = -2.6, -2.4
		double fivePointMid = NumericalFunction.fivePointMidpoint(f, -2.6, h);
		p(fivePointMid);
		p(error(fPrime.apply(-2.6), fivePointMid));
		
		fivePointMid = NumericalFunction.fivePointMidpoint(f, -2.4, h);
		p(fivePointMid);
		p(error(fPrime.apply(-2.4), fivePointMid));
		*/
		
		//5 point endpoint x = -3.0, -2.0
		x0 = -3.0;
		double fivePointEndpoint = NumericalFunction.fivePointEndpoint(f, x0, h);
		p(fivePointEndpoint);
		p(error(fPrime.apply(x0), fivePointEndpoint));
		
		x0 = -2.0;
		fivePointEndpoint = NumericalFunction.fivePointEndpoint(f, x0, -h);
		p(fivePointEndpoint);
		p(error(fPrime.apply(x0), fivePointEndpoint));
		
		/*
		//3 point midpoint x = -2.8, -2.6, -2.4, -2.2
		double threePointMid = 0.0;
		x0 = -2.8;
		threePointMid = NumericalFunction.threePointMidpoint(f, x0, h);
		p(threePointMid);
		p(error(fPrime.apply(x0), threePointMid));
		
		x0 = -2.6;
		threePointMid = NumericalFunction.threePointMidpoint(f, x0, h);
		p(threePointMid);
		p(error(fPrime.apply(x0), threePointMid));
		
		x0 = -2.4;
		threePointMid = NumericalFunction.threePointMidpoint(f, x0, h);
		p(threePointMid);
		p(error(fPrime.apply(x0), threePointMid));
		
		x0 = -2.2;
		threePointMid = NumericalFunction.threePointMidpoint(f, x0, h);
		p(threePointMid);
		p(error(fPrime.apply(x0), threePointMid));
		
		
		//3 point endpoint x= -3.0, -2.0
		double threePointEnd = 0.0;
		x0 = -3.0;
		threePointEnd = NumericalFunction.threePointEndpoint(f, x0, h);
		p(threePointEnd);
		p(error(fPrime.apply(x0), threePointEnd));
		
		x0 = -2.0;
		threePointEnd = NumericalFunction.threePointEndpoint(f, x0, -h);
		p(threePointEnd);
		p(error(fPrime.apply(x0), threePointEnd));
		
		
		double[] x = { -3.0, -2.8, -2.6, -2.4, -2.2, -2.0 };
		for (double xValue : x) {
			System.out.printf("%.11f\n", fPrime.apply(xValue));
		}
		*/
		
	}
	
	private static double error(double expected, double actual) {
		return Math.abs(expected - actual);
	}
	
	private static void p(double x) {
		System.out.printf("%.11f\n", x);
	}
}
