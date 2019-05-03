package tests;

import static org.junit.Assert.*;

import java.util.function.Function;

import org.junit.Test;

import numerical_integration.Integrate;

public class CompositeSimpsonIntegration {

	@Test
	public void test() {
		Function<Double, Double> f = (x) -> {
			return Math.pow(Math.E, x);
		};
		
		double a = 0.0;
		double b = 4.0;
		int n = 4;
		double integral = Integrate.compositeSimpson(a, b, n, f);
		System.out.println(integral);
		assertEquals(53.59819, integral, 1.0);
	}

}
