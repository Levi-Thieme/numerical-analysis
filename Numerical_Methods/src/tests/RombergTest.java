package tests;

import static org.junit.Assert.*;

import java.util.function.Function;

import org.junit.Test;

import numerical_integration.Integrate;

public class RombergTest {

	@Test
	public void testRomberg() {
		Function<Double, Double> f = (x) -> {
			return 1 / (x * Math.log(x));
		};
		double e = Math.E;
		double a = e;
		double b = 2 * e;
		int n = 10;
		double[][] romberg = Integrate.romberg(a, b, f, n, 0.0, 0.0);
		int maxRow = romberg.length - 1;
		int maxColumn = romberg[maxRow].length - 1;
		assertEquals(Math.log(Math.log(2 * e)), romberg[maxRow][maxColumn], 10E-6);
	}
}
