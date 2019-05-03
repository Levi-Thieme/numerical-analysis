package final_exam;

import linear_algebra.LinearSystem;

public class Problem4 {

	public static void main(String[] args) {
		double[][] m = { 
				{5, 1, -6},
				{2, 1, -1},
				{6, 12, 1} };
		double determinant = LinearSystem.determinant3(m);
		System.out.printf("|A| = %f\n", determinant);
		
		double[][] augmentedMatrix = { 
				{5, 1, -6, 7},
				{2, 1, -1, 8},
				{6, 12, 1, 9} };
		double[] x = LinearSystem.gaussPartialPivot(augmentedMatrix);
		double[] expected = { 448.0 / 51.0, -209.0 / 51.0, 93.0 / 17.0 };
		for (int i = 0; i < expected.length; i++) {
			System.out.printf("Error = %.32f\n", Math.abs(x[i] - expected[i]));
			System.out.println(expected[i]);
			System.out.println(x[i]);
		}
	}

}
