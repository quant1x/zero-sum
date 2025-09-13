package main

import (
	"fmt"
	"math"
)

type Point struct {
	X, Y float64
}

// 二次回归函数
func quadraticRegression(points []Point) (func(x float64) float64, [3]float64, error) {
	if len(points) != 3 {
		return nil, [3]float64{}, fmt.Errorf("需要3个点进行二次回归")
	}

	// 构建方程组
	var equations [3][3]float64
	var results [3]float64
	for i, p := range points {
		x, y := p.X, p.Y
		equations[i][0] = x * x
		equations[i][1] = x
		equations[i][2] = 1
		results[i] = y
	}

	// 高斯消元求解
	a, err := solve3x3(equations, results)
	if err != nil {
		return nil, a, err
	}

	// 创建预测函数
	predict := func(x float64) float64 {
		return a[0]*x*x + a[1]*x + a[2]
	}

	return predict, a, nil
}

// 高斯消元法解3x3线性方程组
func solve3x3(A [3][3]float64, b [3]float64) ([3]float64, error) {
	augmented := [3][4]float64{
		{A[0][0], A[0][1], A[0][2], b[0]},
		{A[1][0], A[1][1], A[1][2], b[1]},
		{A[2][0], A[2][1], A[2][2], b[2]},
	}

	for col := 0; col < 3; col++ {
		maxRow := col
		for r := col; r < 3; r++ {
			if math.Abs(augmented[r][col]) > math.Abs(augmented[maxRow][col]) {
				maxRow = r
			}
		}
		augmented[col], augmented[maxRow] = augmented[maxRow], augmented[col]
		pivot := augmented[col][col]
		if math.Abs(pivot) < 1e-9 {
			return [3]float64{}, fmt.Errorf("矩阵奇异")
		}
		for j := col; j < 4; j++ {
			augmented[col][j] /= pivot
		}
		for r := 0; r < 3; r++ {
			if r != col {
				factor := augmented[r][col]
				for j := col; j < 4; j++ {
					augmented[r][j] -= factor * augmented[col][j]
				}
			}
		}
	}

	var solution [3]float64
	for i := 0; i < 3; i++ {
		solution[i] = augmented[i][3]
	}
	return solution, nil
}

// 运动方向分析
func analyzeMotion(a, b, c float64, x float64) (direction string, opening string) {
	// 开口方向判断
	if a > 0 {
		opening = "向上"
	} else {
		opening = "向下"
	}

	// 导数计算
	derivative := 2*a*x + b
	if derivative > 0 {
		direction = "向上"
	} else if derivative < 0 {
		direction = "向下"
	} else {
		direction = "顶点"
	}

	return direction, opening
}

func main() {
	points := []Point{
		{1, 1},
		{2, 2},
		{3, 0},
	}
	points = []Point{
		{X: 0, Y: 6.11},
		{X: 1, Y: 6.18},
		{X: 2, Y: 6.32},
		{X: 3, Y: 6.09},
	}
	predict, coeffs, err := quadraticRegression(points)
	if err != nil {
		fmt.Println("回归失败:", err)
		return
	}
	nextX := len(points)
	for x := range points {
		if nextX < x {
			nextX = x
		}
	}

	// 计算第4个点（假设x=4）
	x4 := float64(nextX + 1)
	y4 := predict(x4)
	fmt.Printf("第四个点坐标: (%.1f, %.2f)\n", x4, y4)

	// 运动方向分析
	direction, opening := analyzeMotion(coeffs[0], coeffs[1], coeffs[2], x4)
	fmt.Printf("曲线开口方向: %s\n", opening)
	fmt.Printf("在x=%.1f处的运动方向: %s\n", x4, direction)
}
