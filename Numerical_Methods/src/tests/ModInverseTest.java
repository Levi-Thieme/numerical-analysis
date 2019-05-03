package tests;

import static org.junit.Assert.*;

import org.junit.Test;

import numerical_functions.NumericalFunction;

public class ModInverseTest {

	@Test
	public void test() {
		long a = 3;
		long m = 11;
		long inverse = NumericalFunction.modInverse(a, m);
		assertEquals(4, inverse);
	}

}
