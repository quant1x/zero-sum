package main

import (
	"fmt"
	"math"
)

type Point struct{ X, Y float64 }

// 目标函数（残差）
func residuals(params []float64, points []Point, res []float64) {
	theta, v0, g := params[0], params[1], params[2]
	cosTheta := math.Cos(theta)
	if cosTheta < 1e-9 {
		cosTheta = 1e-9
	}
	for i, p := range points {
		predictedY := p.X*math.Tan(theta) - (g*p.X*p.X)/(2*v0*v0*cosTheta*cosTheta)
		res[i] = predictedY - p.Y // 注意符号方向
	}
}

// 雅可比矩阵计算
func jacobian(params []float64, points []Point, jac [][]float64) {
	theta, v0, g := params[0], params[1], params[2]
	cosTheta := math.Cos(theta)
	sinTheta := math.Sin(theta)
	cos2Theta := cosTheta * cosTheta
	cos3Theta := cos2Theta * cosTheta
	v0Squared := v0 * v0
	v0Cubed := v0Squared * v0

	for i, p := range points {
		x := p.X
		// df/dθ
		jac[i][0] = (x*(1+math.Tan(theta)*math.Tan(theta)) +
			(g*x*x*sinTheta)/(v0Squared*cos3Theta))

		// df/dv0
		jac[i][1] = -(g * x * x) / (v0Cubed * cos2Theta)

		// df/dg
		jac[i][2] = (x * x) / (2 * v0Squared * cos2Theta)
	}
}

// LM算法实现
func optimizeLM(points []Point, initialParams []float64, maxIter int, tol float64) ([]float64, error) {
	params := make([]float64, len(initialParams))
	copy(params, initialParams)

	lambda := 0.01
	diag := make([]float64, 3)
	res := make([]float64, len(points))
	jac := make([][]float64, len(points))
	for i := range jac {
		jac[i] = make([]float64, 3)
	}

	for iter := 0; iter < maxIter; iter++ {
		residuals(params, points, res)
		jacobian(params, points, jac)

		// 构造正规方程
		JTJ := make([][]float64, 3)
		JTr := make([]float64, 3)
		for i := 0; i < 3; i++ {
			JTJ[i] = make([]float64, 3)
			for j := 0; j < 3; j++ {
				sum := 0.0
				for k := 0; k < len(points); k++ {
					sum += jac[k][i] * jac[k][j]
				}
				JTJ[i][j] = sum
			}
			sum := 0.0
			for k := 0; k < len(points); k++ {
				sum += jac[k][i] * res[k]
			}
			JTr[i] = sum
		}

		// 添加LM阻尼项
		for i := 0; i < 3; i++ {
			diag[i] = JTJ[i][i]
			JTJ[i][i] += lambda
		}

		// 求解增量方程
		delta := make([]float64, 3)
		solve3x3(JTJ, JTr, delta)

		// 尝试更新参数
		newParams := make([]float64, 3)
		copy(newParams, params)
		for i := 0; i < 3; i++ {
			newParams[i] -= delta[i] // 注意符号
		}

		// 计算残差变化
		newRes := make([]float64, len(points))
		residuals(newParams, points, newRes)
		oldNorm := norm(res)
		newNorm := norm(newRes)

		if newNorm < oldNorm {
			copy(params, newParams)
			lambda /= 10
		} else {
			lambda *= 10
		}

		// 收敛检查
		if math.Abs(newNorm-oldNorm) < tol {
			break
		}
	}
	return params, nil
}

func solve3x3(A [][]float64, b []float64, x []float64) {
	// 高斯消去法实现（同前）
	// 此处省略具体实现...
}

func norm(v []float64) float64 {
	sum := 0.0
	for _, val := range v {
		sum += val * val
	}
	return math.Sqrt(sum)
}

func main() {
	points := []Point{
		{1, 1},
		{2, 2},
		{3, 0},
	}
	points = []Point{
		{X: 1, Y: 6.18},
		{X: 2, Y: 6.32},
		{X: 3, Y: 6.09},
	}

	initialParams := []float64{
		math.Atan(1.5), // 初始角度
		5.0,            // 初始速度
		9.8,            // 初始重力
	}

	params, _ := optimizeLM(points, initialParams, 100, 1e-6)
	theta, v0, g := params[0], params[1], params[2]

	fmt.Printf("角度: %.2f°\n", theta*180/math.Pi)
	fmt.Printf("初速度: %.2f m/s\n", v0)
	fmt.Printf("重力加速度: %.2f m/s²\n", g)
}
