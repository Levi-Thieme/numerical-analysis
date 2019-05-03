package tests;

import static org.junit.Assert.*;

import org.junit.Test;

import numerical_functions.NumericalFunction;

public class DividedDifferencesTest {

	@Test
	public void testDividedDiff() {
		double[] x = {1.0, 1.3};
		double[] y = {.7651977, .6200860};
		double dividedDifference = NumericalFunction.dividedDiff(x, y, 0, 1);
		assertEquals(-.4837057, dividedDifference, 1E-7);
	}
	
	@Test
	public void testDividedDifferences() {
		Double[] x = {1.0, 1.3, 1.6, 1.9, 2.2};
		Double[] y = {.7651977, .6200860, .4554022, .2818186, .1103623};
		double[][] diffs = NumericalFunction.dividedDiff(x, y);
		double[] expectedDifferences = { .7651977, -.4837057, -.1087339, .0658784, .0018251 };
		for (int i = 1; i < x.length; i++) {
			System.out.println(diffs[i][i]);
			assertEquals(expectedDifferences[i], diffs[i][i], 1E-7);
		}
	}
	
	@Test
	public void testExtendDividedDifferences() {
		Double[] x = {1.0, 1.3, 1.6, 1.9};
		Double[] y = {.7651977, .6200860, .4554022, .2818186};
		double[][] diffs = NumericalFunction.dividedDiff(x, y);
		double[] expectedDifferences = { .1103623, -.5715210, .0118183, .0680685, .0018251 };
		Double[] newX = {1.0, 1.3, 1.6, 1.9, 2.2};
		Double[] newY = {.7651977, .6200860, .4554022, .2818186, .1103623};
		diffs = NumericalFunction.extendDividedDifferences(newX, newY, diffs);
		for (int i = 1; i < diffs[4].length; i++) {
			System.out.println(diffs[4][i]);
			assertEquals(expectedDifferences[i], diffs[4][i], 1E-7);
		}
	}
}
