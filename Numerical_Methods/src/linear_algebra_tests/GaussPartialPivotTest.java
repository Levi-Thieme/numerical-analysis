package linear_algebra_tests;

import static org.junit.Assert.*;

import org.junit.Test;

import linear_algebra.LinearSystem;

public class GaussPartialPivotTest {
	@Test
	public void test3x3() {
		double expectedX1 = 448.0 / 51.0;
		double expectedX2 = -209.0 / 51.0;
		double expectedX3 = 93.0 / 17.0;
		double[][] augmentedMatrix = { 
				{5, 1, -6, 7},
				{2, 1, -1, 8},
				{6, 12, 1, 9} };
		double[] x = LinearSystem.gaussPartialPivot(augmentedMatrix);
		double accuracy = 1E-4;
		assertEquals(expectedX1, x[0], accuracy);
		assertEquals(expectedX2, x[1], accuracy);
		assertEquals(expectedX3, x[2], accuracy);
		for (double solution : x) {
			System.out.printf("%.18f\n", solution);
		}
	}
	
	
	public void testExample2() {
		double expectedX1 = 10.0;
		double expectedX2 = 1.0;
		double[][] augmentedMatrix = { 
				{ 0.003, 59.14, 59.17},
				{ 5.291, -6.13, 46.78} };
		double[] x = LinearSystem.gaussPartialPivot(augmentedMatrix);
		double accuracy = 1E-4;
		assertEquals(expectedX1, x[0], accuracy);
		assertEquals(expectedX2, x[1], accuracy);
	}

}
