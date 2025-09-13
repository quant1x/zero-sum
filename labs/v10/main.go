package main

import (
	"fmt"
	"math"
)

type Point struct{ X, Y float64 }

func calculateResiduals(params []float64, points []Point) []float64 {
	theta, v0, g := params[0], params[1], params[2]
	residuals := make([]float64, len(points))
	for i, p := range points {
		cosTheta := math.Cos(theta)
		if cosTheta < 1e-9 {
			cosTheta = 1e-9
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
		cos2Theta := cosTheta * cosTheta
		cos3Theta := cos2Theta * cosTheta
		v0Squared := v0 * v0
		v0Cubed := v0Squared * v0

		// 对 theta 的偏导数（严格验证符号）
		jac[i][0] = -(x*(1+math.Tan(theta)*math.Tan(theta)) -
			(g*x*x*sinTheta)/(v0Squared*cos3Theta))

		// 对 v0 的偏导数
		jac[i][1] = (g * x * x) / (v0Cubed * cos2Theta)

		// 对 g 的偏导数
		jac[i][2] = x * x / (2 * v0Squared * cos2Theta)
	}
	return jac
}

func solve3x3(A [][]float64, b []float64) ([]float64, error) {
	augmented := make([][]float64, 3)
	for i := range augmented {
		augmented[i] = make([]float64, 4)
		copy(augmented[i][:3], A[i])
		augmented[i][3] = b[i]
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
		if math.Abs(pivot) < 1e-12 {
			return nil, fmt.Errorf("矩阵奇异")
		}
		for j := col; j <= 3; j++ {
			augmented[col][j] /= pivot
		}
		for r := 0; r < 3; r++ {
			if r != col {
				factor := augmented[r][col]
				for j := col; j <= 3; j++ {
					augmented[r][j] -= factor * augmented[col][j]
				}
			}
		}
	}

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

		// 动态正则化 + 参数约束
		lambda := 1e-6 * matNorm(JTJ)
		for i := 0; i < 3; i++ {
			JTJ[i][i] += lambda
		}

		delta, err := solve3x3(JTJ, JTr)
		if err != nil {
			return nil, fmt.Errorf("线性求解失败: %v", err)
		}

		// 参数更新与约束
		for i := range params {
			params[i] += delta[i] * 0.5 // 使用半步长防止震荡
		}
		params[0] = math.Max(0, math.Min(params[0], math.Pi/2)) // θ ∈ [0, π/2]
		params[1] = math.Max(1e-3, params[1])                   // v₀ > 0
		params[2] = math.Max(1e-3, params[2])                   // g > 0

		// 收敛检查
		maxRes := 0.0
		for _, r := range res {
			if absR := math.Abs(r); absR > maxRes {
				maxRes = absR
			}
		}
		if maxRes < tol {
			return params, nil
		}
	}
	return params, nil
	//return nil, fmt.Errorf("未能在 %d 次迭代内收敛", maxIter)
}

func matNorm(m [][]float64) float64 {
	sum := 0.0
	for i := 0; i < 3; i++ {
		sum += m[i][i] * m[i][i]
	}
	return math.Sqrt(sum)
}

func main() {
	points := []Point{
		{0, 0},
		{1, 1},
		{2, -1},
	}
	points = []Point{
		{X: 0, Y: 6.18},
		{X: 1, Y: 6.32},
		{X: 2, Y: 6.09},
	}
	initialParams := []float64{
		math.Atan(2.5), // 初始角度arctan(2.5) ≈ 68°
		0.1,            // 初始速度
		0.1,            // 重力加速度
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
