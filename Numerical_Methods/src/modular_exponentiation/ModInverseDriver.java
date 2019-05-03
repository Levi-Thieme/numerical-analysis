package modular_exponentiation;

import numerical_functions.NumericalFunction;

public class ModInverseDriver {

	public static void main(String[] args) {
		long a = 3;
		long m = 880;
		long inverse = NumericalFunction.modInverse(a, m);
		System.out.println("Inverse: " + inverse);
	}

}
