package tests;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.function.BiFunction;
import java.util.function.Function;

import org.junit.Test;

import numerical_functions.NumericalFunction;

public class NewtonInterpTest {

	@Test
	public void testNewtonInterp() {
		Double[] x = {1.0, 1.3, 1.6, 1.9, 2.2};
		Double[] y = {.7651977, .6200860, .4554022, .2818186, .1103623};
		ArrayList<BiFunction<Double[], Double[], Double>> coefficientPolynomials = null;
		double[][] dividedDifferences= NumericalFunction.dividedDiff(x, y);
		Function<Double, Double> f = NumericalFunction.newtonInterp(x, y, coefficientPolynomials, dividedDifferences);
		for (int i = 0; i < x.length; i++) {
			double actual = f.apply(x[i]);
			System.out.println(actual);
			assertEquals(y[i], actual, 1E-7);
		}
	}
}
