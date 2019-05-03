package final_exam;


import linear_algebra.LinearSystem;

public class Problem5 {

	public static void main(String[] args) {
		sor();
	}
	
	private static void gaussSeidel() {
		double[][] a = {{3, -1, 1},
				{3, 6, 2},
				{3, 3, 7}};
		double[] x = {0, 0, 0, 0};
		double[] b = {1, 0, 4};
		double[] actualSolutions = { 2.0 / 57.0, -9.0 / 38.0, 25.0 / 38.0};
		double[][] approximatedSolutions = LinearSystem.gaussSeidel(a, b, x, 5);
		
		printMatrix(approximatedSolutions);
		printError(approximatedSolutions, actualSolutions);
	}
	
	private static void jacobi() {
		double[][] a = {{3, -1, 1},
						{3, 6, 2},
						{3, 3, 7}};
		double[] x = {0, 0, 0, 0};
		double[] b = {1, 0, 4};
		double[] actualSolutions = { 2.0 / 57.0, -9.0 / 38.0, 25.0 / 38.0};
		double[][] approximatedSolutions = LinearSystem.jacobi(a, b, x, 5);
		
		printMatrix(approximatedSolutions);
		printError(approximatedSolutions, actualSolutions);
	}
	
	private static void sor() {
		double[][] a = {{3, -1, 1},
				{3, 6, 2},
				{3, 3, 7}};
		double[] x = {0, 0, 0, 0};
		double[] b = {1, 0, 4};
		double[] actualSolutions = { 2.0 / 57.0, -9.0 / 38.0, 25.0 / 38.0};
		double[][] approximatedSolutions = LinearSystem.sor(a, b, x, 5, 1.1, 1E-4);
		printMatrix(approximatedSolutions);
		printError(approximatedSolutions, actualSolutions);
	}
	
	private static void printMatrix(double[][] mat) {
		for (int i = 0; i < mat.length; i++) {
			System.out.printf("X = %.8f Y =  %.8f Z = %.8f\n", mat[i][0], mat[i][1], mat[i][2]);
		}
	}
	
	private static void printError(double[][] mat, double[] actual) {
		for (int i = 0; i < mat.length; i++) {
			System.out.printf("X = %.8f Y = %.8f Z = %.8f\n", Math.abs(mat[i][0] - actual[0]), Math.abs(mat[i][1] - actual[1]), Math.abs(mat[i][2] - actual[2]));
		}
	}

}
