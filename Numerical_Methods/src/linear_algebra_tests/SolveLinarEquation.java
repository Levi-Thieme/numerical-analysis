package linear_algebra_tests;

import static org.junit.Assert.*;

import org.junit.Test;

import linear_algebra.LinearSystem;

public class SolveLinarEquation {

	@Test
	public void test() {
		double[] a = {10, -1, 2, 0};
		double[] x = {0, 0, 0, 0};
		double b = 6;
		int k = 0;
		assertEquals(0.6, LinearSystem.solveForKthVariable(k, a, x, b), 0.0);
	}

}
