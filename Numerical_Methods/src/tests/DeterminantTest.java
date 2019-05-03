package tests;

import static org.junit.Assert.*;

import org.junit.Test;

import linear_algebra.LinearSystem;

public class DeterminantTest {

	@Test
	public void test3x3Determinant() {
		double[][] m = { 
				{6, 1, 1},
				{4, -2, 5},
				{2, 8, 7} };
		double determinant = LinearSystem.determinant3(m);
		assertEquals(-306, determinant, 0.0);
	}

}
