package tests;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.function.BiFunction;
import java.util.function.Function;

import org.junit.Test;

import numerical_functions.NumericalFunction;

public class LagrangeInterpNDegreeTest {
	
	@Test
	public void testLagrangeCoefficientPolynomialsTerm1() {
		Double[] x = {2.0, 2.75, 4.0};
		Double[] xAndK = { 3.0, 0.0};
		ArrayList<BiFunction<Double[], Double[], Double>> L = NumericalFunction.lagrangeCoefficientPolynomials(x.length);
		double approximation = L.get(0).apply(x, xAndK);
		assertEquals(-0.1666666, approximation, .0000001);
	}
	
	@Test
	public void testLagrangeCoefficientPolynomialsTerm2() {
		Double[] x = {2.0, 2.75, 4.0};
		Double[] xAndK = { 3.0, 1.0};
		ArrayList<BiFunction<Double[], Double[], Double>> L = NumericalFunction.lagrangeCoefficientPolynomials(x.length);
		double approximation = L.get(1).apply(x, xAndK);
		assertEquals(1.0666666, approximation, .0000001);
	}
	
	@Test
	public void testLagrangeCoefficientPolynomialsTerm3() {
		Double[] x = {2.0, 2.75, 4.0};
		Double[] xAndK = { 3.0, 2.0};
		ArrayList<BiFunction<Double[], Double[], Double>> L = NumericalFunction.lagrangeCoefficientPolynomials(x.length);
		double approximation = L.get(2).apply(x, xAndK);
		assertEquals(0.1, approximation, .0000001);
	}
	
	@Test
	public void testNDegreeLagrangeInterp() {
		Double[] x = {2.0, 2.75, 4.0};
		double[] y = {.5, .363636, .25};
		ArrayList<BiFunction<Double[], Double[], Double>> L = NumericalFunction.lagrangeCoefficientPolynomials(x.length);
		Function<Double, Double> f = NumericalFunction.lagrangeInterp(x, y, L);
		double approximation = f.apply(3.0);
		System.out.println(approximation);
		assertEquals(0.32955, approximation, .00001);
	}
}
