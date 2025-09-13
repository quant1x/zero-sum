package main

import (
	"fmt"
	"math"
)

type Point struct {
	X, Y float64
}

func calculateResiduals(params []float64, points []Point) []float64 {
	theta, v0, g := params[0], params[1], params[2]
	residuals := make([]float64, len(points))
	for i, p := range points {
		cosTheta := math.Cos(theta)
		if cosTheta == 0 {
			cosTheta = 1e-9 // 防止除以零
		}
		predictedY := p.X*math.Tan(theta) - (g*p.X*p.X)/(2*v0*v0*cosTheta*cosTheta)
		residuals[i] = p.Y - predictedY
	}
	return residuals
}

func calculateJacobian(params []float64, points []Point) [][]float64 {
	theta, v0, g := params[0], params[1], params[2]
	jac := make([][]float64, len(points))
	for i := range jac {
		jac[i] = make([]float64, 3)
	}
	for i, p := range points {
		x := p.X
		cosTheta := math.Cos(theta)
		sinTheta := math.Sin(theta)
		cos3Theta := cosTheta * cosTheta * cosTheta
		v0Squared := v0 * v0

		// 对 theta 的偏导数
		jac[i][0] = (x * (1 + math.Tan(theta)*math.Tan(theta))) -
			(g*x*x*sinTheta)/(v0Squared*cos3Theta)

		// 对 v0 的偏导数
		jac[i][1] = (g * x * x * sinTheta) / (math.Pow(v0, 3) * cosTheta * cosTheta)

		// 对 g 的偏导数
		jac[i][2] = -x * x / (2 * v0Squared * cosTheta * cosTheta)
	}
	return jac
}

func solve3x3(A [][]float64, b []float64) ([]float64, error) {
	// 创建增广矩阵
	augmented := make([][]float64, 3)
	for i := range augmented {
		augmented[i] = make([]float64, 4)
		copy(augmented[i][:3], A[i])
		augmented[i][3] = b[i]
	}

	// 前向消元
	for col := 0; col < 3; col++ {
		// 找主元
		maxRow := col
		for r := col; r < 3; r++ {
			if math.Abs(augmented[r][col]) > math.Abs(augmented[maxRow][col]) {
				maxRow = r
			}
		}
		// 交换行
		augmented[col], augmented[maxRow] = augmented[maxRow], augmented[col]

		pivot := augmented[col][col]
		if math.Abs(pivot) < 1e-12 {
			return nil, fmt.Errorf("矩阵奇异")
		}

		// 归一化当前行
		for j := col; j <= 3; j++ {
			augmented[col][j] /= pivot
		}

		// 消去其他行
		for r := 0; r < 3; r++ {
			if r != col && math.Abs(augmented[r][col]) > 1e-12 {
				factor := augmented[r][col]
				for j := col; j <= 3; j++ {
					augmented[r][j] -= factor * augmented[col][j]
				}
			}
		}
	}

	// 提取解
	x := make([]float64, 3)
	for i := 0; i < 3; i++ {
		x[i] = augmented[i][3]
	}
	return x, nil
}

func solveNewton(params []float64, points []Point, maxIter int, tol float64) ([]float64, error) {
	for iter := 0; iter < maxIter; iter++ {
		res := calculateResiduals(params, points)
		jac := calculateJacobian(params, points)

		// 构造正规方程 J^T*J*Δ = -J^T*res
		JT := make([][]float64, 3)
		for i := 0; i < 3; i++ {
			JT[i] = make([]float64, len(points))
			for j := range points {
				JT[i][j] = jac[j][i]
			}
		}

		JTJ := make([][]float64, 3)
		JTr := make([]float64, 3)
		for i := 0; i < 3; i++ {
			JTJ[i] = make([]float64, 3)
			for j := 0; j < 3; j++ {
				sum := 0.0
				for k := 0; k < len(points); k++ {
					sum += JT[i][k] * jac[k][j]
				}
				JTJ[i][j] = sum
			}
			sum := 0.0
			for k := 0; k < len(points); k++ {
				sum += JT[i][k] * res[k]
			}
			JTr[i] = -sum
		}

		delta, err := solve3x3(JTJ, JTr)
		if err != nil {
			return nil, fmt.Errorf("线性求解失败: %v", err)
		}

		// 更新参数
		for i := range params {
			params[i] += delta[i]
		}

		// 检查收敛
		maxRes := 0.0
		for _, r := range res {
			if math.Abs(r) > maxRes {
				maxRes = math.Abs(r)
			}
		}
		if maxRes < tol {
			return params, nil
		}
	}
	return nil, fmt.Errorf("未能在 %d 次迭代内收敛", maxIter)
}

func main() {
	points := []Point{
		{1, 1},
		{2, 2},
		{3, 0},
	}

	// 初始猜测值（根据物理意义调整）
	initialParams := []float64{
		math.Atan(2.0 / 1.0), // 初始角度根据前两点斜率
		5.0,                  // 初始速度
		9.8,                  // 标准重力加速度
	}

	params, err := solveNewton(initialParams, points, 1000, 1e-8)
	if err != nil {
		fmt.Println("错误:", err)
		return
	}

	theta, v0, g := params[0], params[1], params[2]
	fmt.Printf("角度: %.2f 度\n", theta*180/math.Pi)
	fmt.Printf("初速度: %.2f m/s\n", v0)
	fmt.Printf("重力加速度: %.2f m/s²\n", g)
}
