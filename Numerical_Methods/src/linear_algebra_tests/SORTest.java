package linear_algebra_tests;

import static org.junit.Assert.*;

import org.junit.Test;

import linear_algebra.LinearSystem;

public class SORTest {

	@Test
	public void test() {
		double[][] a = {
				{4, 3, 0}, 
				{3, 4, -1}, 
				{0, -1, 4}};
		double[] b = {24, 30, -24};
		double[] initialSolutions = {1, 1, 1};
		double[][] approximations = LinearSystem.sor(a, b, initialSolutions, 7, 1.25, 1E-4);
		double[] expected = { 3.0, 4.0, -5.0 };
		for (int i = 0; i < approximations[6].length; i++) {
			assertEquals(expected[i], approximations[6][i], 1E-3);
		}
	}

}
