package tests;

import static org.junit.Assert.*;

import org.junit.Test;

import linear_algebra.Norm;

public class MatrixAlgebraTests {

	@Test
	public void test() {
		double[][] m = { {1, 2, -1},
						 {0, 3, -1},
						 {5, -1, 1} };
		
		assertEquals(7, (int)Norm.matrixInfiniteNorm(m));
	}
	
	@Test
	public void test2() {
		double[][] m = { { 2, -1,  0},
						 {-1,  2, -1},
						 { 0,  1,  2} };
		assertEquals(4, (int)Norm.matrixInfiniteNorm(m));
	}

}
