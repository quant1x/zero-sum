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
	jac := make([][]float64, len(points))
	for i := range jac {
		jac[i] = make([]float64, 3)
	}
	for i, p := range points {
		x := p.X
		cosTheta := math.Cos(theta)
		tanTheta := math.Tan(theta)
		secThetaSquared := 1 / (cosTheta * cosTheta)

		// 对 theta 的偏导数（残差函数导数为负）
		dTheta := -(x*(1+tanTheta*tanTheta) - (g*x*x*tanTheta*secThetaSquared)/(v0*v0))
		jac[i][0] = dTheta

		// 对 v0 的偏导数（残差函数导数为负）
		dV0 := -(g * x * x * tanTheta * secThetaSquared) / (math.Pow(v0, 3))
		jac[i][1] = dV0

		// 对 g 的偏导数（残差函数导数为正）
		dG := (x * x) / (2 * v0 * v0 * secThetaSquared)
		jac[i][2] = dG
	}
	return jac
}

func solve3x3(A [][]float64, b []float64) ([]float64, error) {
	a := make([][]float64, 3)
	for i := range a {
		a[i] = make([]float64, 3)
		copy(a[i], A[i])
	}
	x := make([]float64, 3)
	copy(x, b)

	for k := 0; k < 3; k++ {
		maxRow := k
		max := math.Abs(a[k][k])
		for i := k + 1; i < 3; i++ {
			if abs := math.Abs(a[i][k]); abs > max {
				max = abs
				maxRow = i
			}
		}
		a[k], a[maxRow] = a[maxRow], a[k]
		x[k], x[maxRow] = x[maxRow], x[k]

		for i := k + 1; i < 3; i++ {
			factor := a[i][k] / a[k][k]
			x[i] -= factor * x[k]
			for j := k; j < 3; j++ {
				a[i][j] -= factor * a[k][j]
			}
		}
	}

	for i := 2; i >= 0; i-- {
		x[i] /= a[i][i]
		for j := i - 1; j >= 0; j-- {
			x[j] -= a[j][i] * x[i]
		}
	}

	return x, nil
}

func solveNewton(params []float64, points []Point, maxIter int, tol float64) ([]float64, error) {
	for i := 0; i < maxIter; i++ {
		res := calculateResiduals(params, points)
		jac := calculateJacobian(params, points)

		H := make([][]float64, 3)
		for i := range H {
			H[i] = make([]float64, 3)
		}
		b := make([]float64, 3)

		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				sum := 0.0
				for k := 0; k < len(points); k++ {
					sum += jac[k][i] * jac[k][j]
				}
				H[i][j] = sum
			}
			b[i] = 0.0
			for k := 0; k < len(points); k++ {
				b[i] -= jac[k][i] * res[k]
			}
		}

		delta, err := solve3x3(H, b)
		if err != nil {
			return nil, err
		}

		for j := range params {
			params[j] += delta[j]
		}

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

	// 调整初始猜测值（更接近真实解）
	initialParams := []float64{math.Atan(4), 9.0, 9.8}

	params, err := solveNewton(initialParams, points, 1000, 1e-8)
	if err != nil {
		fmt.Println("错误:", err)
		return
	}

	theta, v0, g := params[0], params[1], params[2]
	fmt.Printf("角度: %.2f 度\n", theta*180/math.Pi)
	fmt.Printf("初速度: %.2f\n", v0)
	fmt.Printf("重力加速度: %.2f\n", g)
}
