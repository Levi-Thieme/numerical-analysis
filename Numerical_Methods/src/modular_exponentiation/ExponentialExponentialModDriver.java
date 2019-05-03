package modular_exponentiation;

import numerical_functions.NumericalFunction;

public class ExponentialExponentialModDriver {

	public static void main(String[] args) {
		double x = 4;
		int y = 2;
		int z = 2006;
		int p = 31;
		double result = NumericalFunction.exponentialExponentialModulus(x, y, z, p);
		System.out.println("Result: " + result + "\n");
	}
}