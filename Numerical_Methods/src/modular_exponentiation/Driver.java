package modular_exponentiation;

import numerical_functions.NumericalFunction;

public class Driver {

	public static void main(String[] args) {
		double base = 4;
		int pow =(int) Math.pow(2, 26);
		int mod = 31;
		double result = NumericalFunction.exponentialModulus(base, pow, mod);
		System.out.println("Result: " + result + "\n");
		System.out.println(NumericalFunction.exponentialModulus(741, 587, 943));
	}
}
