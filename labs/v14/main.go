package main

import (
	"fmt"
	"math"
)

type Point struct {
	X, Y float64
}

// 二次回归函数（支持≥3个点）
func quadraticRegression(points []Point) (func(x float64) float64, [3]float64, error) {
	n := len(points)
	if n < 3 {
		return nil, [3]float64{}, fmt.Errorf("需要至少3个点")
	}

	// 构建设计矩阵和目标向量
	X := make([][]float64, n)
	y := make([]float64, n)
	for i, p := range points {
		X[i] = []float64{p.X * p.X, p.X, 1}
		y[i] = p.Y
	}

	// 计算正规方程 X^T·X·β = X^T·y
	XT := transpose(X)
	XTX := multiplyMatrices(XT, X)
	XTy := multiplyMatrixVector(XT, y)

	// 解正规方程
	coeffs, err := solveGauss(XTX, XTy)
	if err != nil {
		return nil, [3]float64{}, err
	}

	// 创建预测函数
	predict := func(x float64) float64 {
		return coeffs[0]*x*x + coeffs[1]*x + coeffs[2]
	}

	return predict, [3]float64{coeffs[0], coeffs[1], coeffs[2]}, nil
}

// 矩阵转置
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

// 矩阵乘法
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

// 矩阵向量乘法
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

// 高斯消元法解线性方程组
func solveGauss(A [][]float64, b []float64) ([]float64, error) {
	n := len(A)
	augmented := make([][]float64, n)
	for i := range augmented {
		augmented[i] = append(A[i], b[i])
	}

	for col := 0; col < n; col++ {
		// 选主元
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

		// 归一化
		for j := col; j <= n; j++ {
			augmented[col][j] /= pivot
		}

		// 消元
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

// 运动分析
func analyzeMotion(a, b, c float64, x float64) (direction string, opening string) {
	// 开口方向判断
	if a > 0 {
		opening = "向上"
	} else {
		opening = "向下"
	}

	// 导数判断运动方向
	derivative := 2*a*x + b
	switch {
	case derivative > 1e-6:
		direction = "向上"
	case derivative < -1e-6:
		direction = "向下"
	default:
		direction = "顶点"
	}
	return
}

func main() {
	// 输入点（可扩展）
	points := []Point{
		{1, 1},
		{2, 2},
		{3, 0},
		{4, -5}, // 新增第四个点验证
	}
	points = []Point{
		{X: 0, Y: 5.64},
		{X: 1, Y: 6.18},
		{X: 2, Y: 6.32},
		{X: 3, Y: 6.09},
	}
	// 执行二次回归
	predict, coeffs, err := quadraticRegression(points)
	if err != nil {
		fmt.Println("回归失败:", err)
		return
	}
	fmt.Printf("回归方程: y = %.2fx² + %.2fx + %.2f\n", coeffs[0], coeffs[1], coeffs[2])
	nextX := 1
	for x := range points {
		if nextX < x {
			nextX = x
		}
	}
	// 预测第四个点（x=5）
	xNew := float64(nextX + 1)
	yNew := predict(xNew)
	fmt.Printf("预测点(%.1f): y = %.2f\n", xNew, yNew)

	// 运动分析
	direction, opening := analyzeMotion(coeffs[0], coeffs[1], coeffs[2], xNew)
	fmt.Printf("曲线开口方向: %s\n", opening)
	fmt.Printf("在x=%.1f处的运动方向: %s\n", xNew, direction)
}
