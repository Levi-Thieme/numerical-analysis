package tests;

import static org.junit.Assert.*;

import org.junit.Test;

import numerical_functions.NumericalFunction;

public class ExtendedGCDTest {

	@Test
	public void testExtendedGCD() {
		long a = 35;
		long b = 15;
		long[] t = {1, 1};
		long gcd = NumericalFunction.gcdExtended(a, b, t);
		assertEquals(1, t[0]);
		assertEquals(-2, t[1]);
		assertEquals(gcd, 5);
	}

}
