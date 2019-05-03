package linear_algebra;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.BiFunction;
import java.util.function.Function;

public class LinearSystem {
	
	
	/**
	 * Uses Gaussian elimination and backward substitution with partial pivoting
	 * to solve a system of linear equations.
	 * 
	 * @param m an augmented matrix
	 * @return An array containing the solutions of the linear system.
	 */
	public static double[] gaussPartialPivot(double[][] m) {
		int maxRow = m.length - 1;
		int maxCol = m[0].length - 1;
		double[] x = new double[m.length];
		
		//elimination process
		//for each column, use length - 1 for augmented matrix
		for (int i = 0; i < m.length - 1; i++) {
			//pivot current row to the maximum of rows i...n
			double[] current = m[i];
			int pivotIndex = i;
			//find the maximum rows below current and pivot
			for (int j = i + 1; j < m.length; j++) {
				if (Math.abs(m[j][i]) > Math.abs(current[i])) {
					pivotIndex = j;
				}
			}
			if (m[pivotIndex][i] == 0) {
				System.out.println("No unique solution exists!");
				return null;
			} 
			else if (m[pivotIndex] != current) {
				//swap the rows, current index is i
				double[] temp = m[i];
				m[i] = m[pivotIndex];
				m[pivotIndex] = temp;
				current = m[i];
				System.out.printf("Row %d pivots with row %d\n", i + 1, pivotIndex + 1);
			}
			//elimination step
			for (int k = i + 1; k < m.length; k++) {
				final double multiplier = m[k][i] / current[i];
				BiFunction<Double, Double, Double> op = (y, z)-> { return y - (multiplier * z); };
				m[k] = rowOperations(m[k], m[i], op);
			}
		}
		if (m[m.length - 1][m[m.length - 1].length - 1] == 0) {
			System.out.println("No unique solution exists!");
			return null;
		}
		//backward substitution
		x[maxRow] = m[maxRow][maxCol] / m[maxRow][maxCol - 1];
		
		for (int i = maxRow - 1; i >= 0; i--) {
			double sum = 0.0;
			for (int j = i + 1; j <= maxRow; j++) {
				sum += m[i][j] * x[j];
			}
			//Xi = Bi - sum / m[i][i];
			x[i] = (m[i][maxCol] - sum) / m[i][i];
		}
		return x;
	}
	
	/**
	 * Performs a row operation specified by the BiFunction op.
	 * @param x double array
	 * @param y double array
	 * @param op A BiFunction accepting two arguments that performs a binary operation on the arguments and returns the result.
	 * @return A new double array containing the results of op performed on x and y.
	 */
	public static double[] rowOperations(double[] x, double[] y, BiFunction<Double, Double, Double> op) {
		int length = Math.min(x.length, y.length);
		double[] transformedRow = new double[Math.max(x.length, y.length)];
		for (int i = 0; i < length; i++) {
			transformedRow[i] = op.apply(x[i], y[i]);
		}
		return transformedRow;
	}
	
	
	/**
	 * Computes the determinant of a 2x2 matrix.
	 * Precondition: m must be a 2x2 matrix.
	 * @param m 2x2 matrix
	 * @return the determinant of m
	 */
	public static double determinant2(double[][] m) {
		// |m| = ad - bc
		double a = m[0][0];
		double b = m[0][1];
		double c = m[1][0];
		double d = m[1][1];
		return a * d - b * c;
	}
	
	/**
	 * Computes the determinant of a 3x3 matrix.
	 * Precondition: m must be a 3x3 matrix.
	 * @param m 3x3 matrix
	 * @return the determinant of m
	 */
	public static double determinant3(double[][] m) {
		/*
		 *  |m| = a|a'| - b|b'| + c|c'|, where k' is formed by removing all entries in the 
		 *  same row and column as k and using the remaining entries to form a new matrix.
		 */
		double a = m[0][0];
		double[][] aPrimeMat = { {m[1][1], m[1][2]}, 
								{ m[2][1], m[2][2] } };
		double aPrime = determinant2(aPrimeMat);
		
		double b = m[0][1];
		double[][] bPrimeMat = { {m[1][0], m[1][2] },
								{ m[2][0], m[2][2] } };
		double bPrime = determinant2(bPrimeMat);
		
		double c = m[0][2];
		double[][] cPrimeMat = { {m[1][0], m[1][1] },
								{ m[2][0], m[2][1] } };
		double cPrime = determinant2(cPrimeMat);
		
		return a * aPrime - b * bPrime + c * cPrime;
	}

	/**
	 * Uses the Jacobi iterative technique to solve a linear system of equations.
	 * 
	 * @param a double matrix
	 * @param b double array of solutions
	 * @param initialSolutions double array of initial solutions
	 * @param tolerance The desired accuracy of the approximated solutions to the actual solutions.
	 * @return A double array containing the approximated solutions.
	 */
	public static double[] jacobi(double[][] a, double[] b, double[] initialSolutions, double tolerance) {
		double[] previousSolutions = initialSolutions;
		double previousMaxNorm = Norm.infiniteNorm(initialSolutions);
		double[] currentSolutions = new double[b.length];
		//solve each row for its diagonal variable
		boolean done = false;
		while (!done) {
			for (int k = 0; k < a.length; k++) {
				//solve a[k] for x[k]
				currentSolutions[k] = solveForKthVariable(k, a[k], previousSolutions, b[k]);				
			}
			//test for the difference in the infinite norms of the previous and 
			//current solution arrays being less than tolerance
			double currentMaxNorm = Norm.infiniteNorm(currentSolutions);
			if (Math.abs(currentMaxNorm - previousMaxNorm) < tolerance) {
				done = true;
			}
			previousSolutions = currentSolutions;
			previousMaxNorm = currentMaxNorm;
		}
		return currentSolutions;
	}
	
	/**
	 * Uses the Jacobi iterative technique to solve a linear system of equations.
	 * 
	 * @param a double matrix
	 * @param b double array of solutions
	 * @param initialSolutions double array of initial solutions
	 * @param maxIterations The maximum number of iterations to perform.
	 * @return A double matrix containing the approximated solutions.
	 */
	public static double[][] jacobi(double[][] a, double[] b, double[] initialSolutions, int maxIterations) {
		double[] previousSolutions = initialSolutions;
		double[][] solutions = new double[maxIterations][b.length];
		for (int i = 0; i < maxIterations; i++) {
			for (int k = 0; k < solutions[k].length; k++) {
				//solve a[k] for x[k]
				solutions[i][k] = solveForKthVariable(k, a[k], previousSolutions, b[k]);			
			}
			previousSolutions = solutions[i];
		}
		return solutions;
	}
	
	/**
	 * Solves a linear equation for the kth variable.
	 * 
	 * @param k The index specifying which variable for which to solve.
	 * @param a Coefficients of the equation.
	 * @param x Variables of the equation.
	 * @param b Solution of the equation.
	 * @return The value of the kth variable in terms of the other coefficients, variables, and solution.
	 */
	public static double solveForKthVariable(int k, double[] a, double[] x, double b) {
		double total = 0.0;
		for (int j = 0; j < a.length; j++) {
			if (j != k) {
				total += (a[j] * x[j]);
			}
		}
		return (b - total) / a[k];
	}
	
	/**
	 * Calculates the sum of the product at each entry for [start...end]
	 * @param a
	 * @param x
	 * @param start
	 * @param end
	 * @return 
	 */
	private static double productSum(double[] a, double[] x, int start, int end) {
		double total = 0.0;
		for (int i = start; i < end; i++) {
			total += a[i] * x[i];
		}
		return total;
	}
	
	/**
	 * Performs the Gauss-Seidel calculation for x[i]
	 * @param i
	 * @param a
	 * @param x
	 * @param xPrevious
	 * @param b
	 * @return
	 */
	public static double gaussSeidelVarSolution(int i, double[][] a, double[] x, double[] xPrevious, double[] b) {
		double currentTotal = productSum(a[i], x, 0, i);
		double previousTotal = productSum(a[i], xPrevious, i + 1, a[i].length);
		return (b[i] - currentTotal - previousTotal) / a[i][i];
	}
	
	/**
	 * Gauss-Seidel method for solving a system of linear equations.
	 * 
	 * @param a coefficient matrix
	 * @param b solution vector
	 * @param initialSolutions initial solution vector
	 * @param maxIterations The maximum number of iterations to perform
	 * @return A matrix containing the solution approximations at each iteration
	 */
	public static double[][] gaussSeidel(double[][] a, double[] b, double[] initialSolutions, int maxIterations) {
		double[][] x = new double[maxIterations][b.length];
		//System.out.println(initialSolutions);
		//solve each row for its diagonal variable
		for (int k = 0; k < maxIterations; k++) {
			for (int i = 0; i < a.length; i++) {
				x[k][i] = gaussSeidelVarSolution(i, a, x[k], initialSolutions, b);
			}	
			initialSolutions = Arrays.copyOf(x[k], initialSolutions.length);
			//System.out.println(initialSolutions);
		}
		return x;
	}
	
	/**
	 * Successive Over Relaxation method for solving systems of linear equations.
	 * 
	 * @param a Coefficient matrix
	 * @param b Solution vector
	 * @param initialSolutions Vector of initial solutions
	 * @param N Maximum number of iterations
	 * @param w Weight for over relaxation. 1 < w < 2
	 * @param tolerance The desired difference between successive solution approximation norms.
	 * @return A matrix containing approximation vectors for iterations 0...N
	 */
	public static double[][] sor(double[][] a, double[] b, double[] initialSolutions, int N, double w, double tolerance) {
		int n = b.length;
		double[][] x = new double[N][n];
		int k = 0;
		while (k < N) {
			for (int i = 0; i < n; i++) {
				x[k][i] = (1.0 - w) * initialSolutions[i] + 
						(1.0 / a[i][i]) * (w * ( - productSum(a[i], x[k], 0, i) - productSum(a[i], initialSolutions, i + 1, n) + b[i]));
			}
			if (Norm.euclideanNorm(x[k]) - Norm.euclideanNorm(initialSolutions) < tolerance) {
				System.out.printf("Tolerance satisfied at iteration %d\n", k + 1);
			}
			for (int i = 0; i < n; i++) {
				initialSolutions[i] = x[k][i];
			}
			k += 1;
		}
		return x;
	}
	
	/**
	 * Sums the entries from start to end inclusively.
	 * @param values
	 * @param start
	 * @param end
	 * @return The sum of values[start...end]
	 */
	public static double sum(double[] values, int start, int end) {
		double total = 0.0;
		for (int i = start; i <= end; i++) {
			total += values[i];
		}
		return total;
	}
	
}






































