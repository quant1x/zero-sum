// Package regression 提供二次回归分析，用于将数据点拟合到二次多项式 y = ax² + bx + c 并分析曲线属性。
package regression

import (
	"errors"
	"math"
)

// Point 表示具有 X 和 Y 坐标的二维点。
type Point struct {
	X, Y float64
}

// QuadraticRegression 对一组点执行二次回归。
// 需要至少 3 个点，返回预测函数、系数 [a, b, c] 和任何错误。
func QuadraticRegression(points []Point) (func(x float64) float64, [3]float64, error) {
	n := len(points)
	if n < 3 {
		return nil, [3]float64{}, errors.New("需要至少3个点")
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

// transpose 矩阵转置
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

// multiplyMatrices 矩阵乘法
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

// multiplyMatrixVector 矩阵向量乘法
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

// solveGauss 高斯消元法解线性方程组
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
			return nil, errors.New("矩阵奇异")
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

// AnalyzeMotion 分析由系数 a, b, c 定义的二次曲线在点 x 处的运动方向和开口方向。
// 返回开口方向（"向上" 或 "向下"）和运动方向（"向上"、"向下" 或 "顶点"）。
func AnalyzeMotion(a, b, c float64, x float64) (direction string, opening string) {
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