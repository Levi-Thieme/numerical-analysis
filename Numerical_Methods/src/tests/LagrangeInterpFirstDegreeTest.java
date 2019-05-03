package tests;

import static org.junit.Assert.*;

import org.junit.Test;

import numerical_functions.NumericalFunction;

public class LagrangeInterpFirstDegreeTest {

	@Test
	public void testLagrangeInterpX0() {
		double[] xVals = {2, 4};
		double[] yVals = {5, 1};
		double x = xVals[0];
		double result = NumericalFunction.lagrangeInterpFirstDegree(x, xVals, yVals);
		System.out.println(result);
		assertEquals(yVals[0], result, 0.0000001);
	}
	
	@Test
	public void testLagrangeInterpX1() {
		double[] xVals = {2, 4};
		double[] yVals = {5, 1};
		double x = xVals[1];
		double result = NumericalFunction.lagrangeInterpFirstDegree(x, xVals, yVals);
		System.out.println(result);
		assertEquals(yVals[1], result, 0.0000001);
	}
}
