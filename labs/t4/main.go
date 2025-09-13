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
		predictedY := p.X*math.Tan(theta) - (g*p.X*p.X)/(2*v0*v0*math.Cos(theta)*math.Cos(theta))
		residuals[i] = p.Y - predictedY
	}
	return residuals
}

func calculateJacobian(params []float64, points []Point) [][]float64 {
	theta, v0, g := params[0], params[1], params[2]
	jacobian := make([][]float64, len(points))
	for i := range jacobian {
		jacobian[i] = make([]float64, 3)
	}
	for i, p := range points {
		x := p.X
		cosTheta := math.Cos(theta)
		tanTheta := math.Tan(theta)
		secThetaSquared := 1 / (cosTheta * cosTheta)

		// 对 theta 的偏导数
		dTheta := x*(1+tanTheta*tanTheta) - (g*x*x*2*tanTheta*secThetaSquared)/(2*v0*v0)
		jacobian[i][0] = dTheta

		// 对 v0 的偏导数
		dV0 := (g * x * x) / (math.Pow(v0, 3) * secThetaSquared)
		jacobian[i][1] = dV0

		// 对 g 的偏导数
		dG := -x * x / (2 * v0 * v0 * secThetaSquared)
		jacobian[i][2] = dG
	}
	return jacobian
}

func solveNewton(params []float64, points []Point, maxIter int, tol float64) ([]float64, error) {
	for i := 0; i < maxIter; i++ {
		res := calculateResiduals(params, points)
		jac := calculateJacobian(params, points)

		// 使用伪逆法求解线性方程组 J * delta = -res
		delta := make([]float64, 3)
		for j := 0; j < 3; j++ {
			sum := 0.0
			for k := 0; k < 3; k++ {
				sum += jac[0][k]*jac[0][j] + jac[1][k]*jac[1][j] + jac[2][k]*jac[2][j]
			}
			delta[j] = -res[0]*jac[0][j] - res[1]*jac[1][j] - res[2]*jac[2][j]
			delta[j] /= sum
		}

		// 更新参数
		for j := range params {
			params[j] += delta[j]
		}

		// 检查收敛性
		maxResidual := 0.0
		for _, r := range res {
			absR := math.Abs(r)
			if absR > maxResidual {
				maxResidual = absR
			}
		}
		if maxResidual < tol {
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

	// 初始猜测值：theta=45度，v0=5，g=9.8
	initialParams := []float64{math.Pi / 4, 5.0, 9.8}
	maxIter := 1_000_0000
	params, err := solveNewton(initialParams, points, maxIter, 1e-6)
	if err != nil {
		fmt.Println("错误:", err)
		return
	}

	theta, v0, g := params[0], params[1], params[2]
	fmt.Printf("角度: %.2f 度\n", theta*180/math.Pi)
	fmt.Printf("初速度: %.2f\n", v0)
	fmt.Printf("重力加速度: %.2f\n", g)
}
