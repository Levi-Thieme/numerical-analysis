package tests;

import static org.junit.Assert.*;

import org.junit.Test;

import numerical_functions.NumericalFunction;

public class GreatestCommonDivisorTest {
	
	@Test
	public void testFirstZero() {
		long a = 0;
		long b = 192;
		long greatestDivisor = NumericalFunction.gcd(a, b);
		assertEquals(greatestDivisor, b);
	}
	
	@Test
	public void testSecondZero() {
		long a = 270;
		long b = 0;
		long greatestDivisor = NumericalFunction.gcd(a, b);
		assertEquals(greatestDivisor, a);
	}
	
	@Test
	public void gcdTestSmall() {
		long a = 270;
		long b = 192;
		long greatestDivisor = NumericalFunction.gcd(a, b);
		assertEquals(greatestDivisor, 6);
	}
	
	@Test
	public void gcdTestLarge() {
		long a = 7966496;
		long b = 314080416;
		long greatestDivisor = NumericalFunction.gcd(a, b);
		assertEquals(greatestDivisor, 32);
	}
}
