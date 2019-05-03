package linear_algebra_tests;

import static org.junit.Assert.*;

import org.junit.Test;

import linear_algebra.LinearSystem;

public class JacobiTest {

	@Test
	public void test() {
		double[][] a = {  {10, -1, 2, 0},
						{ -1, 11, -1, 3 },
						{ 2, -1, 10, -1 },
						{ 0, 3, -1, 8} };
		double[] x = {0, 0, 0, 0};
		double[] b = {6, 25, -11, 15};
		double[] expectedSolutions = {1, 2, -1, 1};
		double[][] actualSolutions = LinearSystem.jacobi(a, b, x, 5);
		for (double solution : actualSolutions[actualSolutions.length - 1]) {
			System.out.printf("%.8f\n", solution);
		}
		for (int i = 0; i < actualSolutions.length; i++) {
			assertEquals(expectedSolutions[i], actualSolutions[actualSolutions.length - 1][i], 1E-3);
		}
	}

}
