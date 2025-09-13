package main

import (
	"fmt"
	"math"
)

type Point struct {
	X, Y float64
}

func polynomialRegression(points []Point, degree int) ([]float64, error) {
	n := len(points)
	if n < degree+1 {
		return nil, fmt.Errorf("需要至少 %d 个点", degree+1)
	}

	// 构建设计矩阵和目标向量
	X := make([][]float64, n)
	y := make([]float64, n)
	for i, p := range points {
		X[i] = make([]float64, degree+1)
		for j := range X[i] {
			X[i][j] = math.Pow(p.X, float64(j))
		}
		y[i] = p.Y
	}

	// 计算正规方程
	XT := transpose(X)
	XTX := multiplyMatrices(XT, X)
	XTy := multiplyMatrixVector(XT, y)

	// 解线性方程组
	coeffs, err := solveGauss(XTX, XTy)
	if err != nil {
		return nil, err
	}

	return coeffs, nil
}

func transpose(matrix [][]float64) [][]float64 {
	rows := len(matrix)
	cols := len(matrix[0])
	transposed := make([][]float64, cols)
	for i := range transposed {
		transposed[i] = make([]float64, rows)
		for j := range transposed[i] {
			transposed[i][j] = matrix[j][i]
		}
	}
	return transposed
}

func multiplyMatrices(a, b [][]float64) [][]float64 {
	rowsA := len(a)
	colsA := len(a[0])
	colsB := len(b[0])
	result := make([][]float64, rowsA)
	for i := range result {
		result[i] = make([]float64, colsB)
		for k := 0; k < colsA; k++ {
			for j := 0; j < colsB; j++ {
				result[i][j] += a[i][k] * b[k][j]
			}
		}
	}
	return result
}

func multiplyMatrixVector(a [][]float64, v []float64) []float64 {
	rowsA := len(a)
	colsA := len(a[0])
	result := make([]float64, rowsA)
	for i := 0; i < rowsA; i++ {
		for j := 0; j < colsA; j++ {
			result[i] += a[i][j] * v[j]
		}
	}
	return result
}

func solveGauss(A [][]float64, b []float64) ([]float64, error) {
	n := len(A)
	augmented := make([][]float64, n)
	for i := range augmented {
		augmented[i] = append(A[i], b[i])
	}

	for col := 0; col < n; col++ {
		maxRow := col
		for r := col; r < n; r++ {
			if math.Abs(augmented[r][col]) > math.Abs(augmented[maxRow][col]) {
				maxRow = r
			}
		}
		augmented[col], augmented[maxRow] = augmented[maxRow], augmented[col]
		pivot := augmented[col][col]
		if math.Abs(pivot) < 1e-9 {
			return nil, fmt.Errorf("矩阵奇异")
		}
		for j := col; j <= n; j++ {
			augmented[col][j] /= pivot
		}
		for r := 0; r < n; r++ {
			if r != col {
				factor := augmented[r][col]
				for j := col; j <= n; j++ {
					augmented[r][j] -= factor * augmented[col][j]
				}
			}
		}
	}

	x := make([]float64, n)
	for i := range x {
		x[i] = augmented[i][n]
	}
	return x, nil
}

func main() {
	points := []Point{
		{1, 1},
		{2, 2},
		{3, 0},
	}

	coeffs, err := polynomialRegression(points, 2)
	if err != nil {
		fmt.Println("错误:", err)
		return
	}

	fmt.Printf("二次回归方程: y = %.2fx² + %.2fx + %.2f\n",
		coeffs[2], coeffs[1], coeffs[0])
}
