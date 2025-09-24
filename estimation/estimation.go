// Package estimation 提供参数估计算法，使用Levenberg-Marquardt方法反推抛物线物理参数。
package estimation

import (
	"errors"
	"math"
)

// Point 表示二维点。
type Point struct {
	X, Y float64
}

// ProjectileSolver 抛物线参数求解器。
type ProjectileSolver struct {
	Points  [3]Point // 三个观测点
	Tol     float64  // 收敛容差
	MaxIter int      // 最大迭代次数
}

// PhysicsParams 物理参数。
type PhysicsParams struct {
	V0    float64 // 初速度
	Theta float64 // 抛射角（弧度）
	G     float64 // 重力加速度
}

// Solve 执行LM算法求解参数。
func (s *ProjectileSolver) Solve(initial PhysicsParams) (PhysicsParams, error) {
	p := initial
	lambda := 1e-3 // 阻尼系数
	prevError := math.MaxFloat64

	for iter := 0; iter < s.MaxIter; iter++ {
		J, F := s.jacobian(p)
		JTJ := multiplyTransposed(J, J)
		JTF := multiplyTransposedVec(J, F)

		// 添加阻尼项：JTJ + λ*I
		for i := 0; i < 3; i++ {
			JTJ[i][i] *= (1 + lambda)
		}

		delta, ok := solveLinearSystem(JTJ, JTF)
		if !ok {
			return p, errors.New("矩阵奇异")
		}

		// 临时更新参数
		tempParams := PhysicsParams{
			V0:    p.V0 - delta[0],
			Theta: p.Theta - delta[1],
			G:     p.G - delta[2],
		}

		// 计算新误差
		newF := s.residuals(tempParams)
		newError := norm(newF)

		// 误差下降则接受更新，否则增大阻尼
		if newError < prevError {
			p = tempParams
			prevError = newError
			lambda *= 0.1
		} else {
			lambda *= 10
		}

		// 物理约束
		p.V0 = math.Max(p.V0, 0.1)
		p.G = math.Max(p.G, 0.1)
		p.Theta = math.Mod(math.Abs(p.Theta), math.Pi/2)

		// 收敛检查
		if newError < s.Tol {
			return p, nil
		}
	}
	return p, errors.New("超过最大迭代次数")
}

// residuals 计算残差向量。
func (s *ProjectileSolver) residuals(p PhysicsParams) [3]float64 {
	var F [3]float64
	for i, pt := range s.Points {
		x, y := pt.X, pt.Y
		F[i] = y - (x*math.Tan(p.Theta) - (p.G*x*x)/(2*p.V0*p.V0*math.Pow(math.Cos(p.Theta), 2)))
	}
	return F
}

// jacobian 计算雅可比矩阵（数值微分法）。
func (s *ProjectileSolver) jacobian(p PhysicsParams) ([3][3]float64, [3]float64) {
	epsilon := 1e-6
	F := s.residuals(p)
	J := [3][3]float64{}

	// 对每个参数扰动计算偏导
	for i := 0; i < 3; i++ {
		perturbed := p
		switch i {
		case 0:
			perturbed.V0 += epsilon
		case 1:
			perturbed.Theta += epsilon
		case 2:
			perturbed.G += epsilon
		}

		F_perturbed := s.residuals(perturbed)
		for j := 0; j < 3; j++ {
			J[j][i] = (F_perturbed[j] - F[j]) / epsilon
		}
	}
	return J, F
}

// multiplyTransposed 计算A^T * B。
func multiplyTransposed(A [3][3]float64, B [3][3]float64) [3][3]float64 {
	var result [3][3]float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			sum := 0.0
			for k := 0; k < 3; k++ {
				sum += A[k][i] * B[k][j]
			}
			result[i][j] = sum
		}
	}
	return result
}

// multiplyTransposedVec 计算A^T * v。
func multiplyTransposedVec(A [3][3]float64, v [3]float64) [3]float64 {
	var result [3]float64
	for i := 0; i < 3; i++ {
		sum := 0.0
		for j := 0; j < 3; j++ {
			sum += A[j][i] * v[j]
		}
		result[i] = sum
	}
	return result
}

// solveLinearSystem 解线性方程组（伴随矩阵法）。
func solveLinearSystem(J [3][3]float64, F [3]float64) ([3]float64, bool) {
	// 计算行列式
	det := J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1]) -
		J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0]) +
		J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0])

	if math.Abs(det) < 1e-20 { // 防止除以零
		return [3]float64{}, false
	}

	// 计算伴随矩阵
	adjugate := [3][3]float64{
		{ // 第一行
			(J[1][1]*J[2][2] - J[1][2]*J[2][1]),
			-(J[0][1]*J[2][2] - J[0][2]*J[2][1]),
			(J[0][1]*J[1][2] - J[0][2]*J[1][1]),
		},
		{ // 第二行
			-(J[1][0]*J[2][2] - J[1][2]*J[2][0]),
			(J[0][0]*J[2][2] - J[0][2]*J[2][0]),
			-(J[0][0]*J[1][2] - J[0][2]*J[1][0]),
		},
		{ // 第三行
			(J[1][0]*J[2][1] - J[1][1]*J[2][0]),
			-(J[0][0]*J[2][1] - J[0][1]*J[2][0]),
			(J[0][0]*J[1][1] - J[0][1]*J[1][0]),
		},
	}

	// 计算逆矩阵 = 伴随矩阵 / 行列式
	inv := [3][3]float64{}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			inv[i][j] = adjugate[j][i] / det // 转置伴随矩阵得到逆矩阵
		}
	}

	// 计算增量 delta = -inv * F
	delta := [3]float64{
		-(inv[0][0]*F[0] + inv[0][1]*F[1] + inv[0][2]*F[2]),
		-(inv[1][0]*F[0] + inv[1][1]*F[1] + inv[1][2]*F[2]),
		-(inv[2][0]*F[0] + inv[2][1]*F[1] + inv[2][2]*F[2]),
	}

	return delta, true
}

// norm 计算向量范数。
func norm(v [3]float64) float64 {
	return math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}

// EstimateInitial 自动估算初始值。
func EstimateInitial(points [3]Point) PhysicsParams {
	// 使用更稳定的初始值估算
	// 基于抛物线公式估算参数

	p2 := points[1] // 使用中间点

	// 简单估算：假设θ=45°，g=9.8
	theta_est := math.Pi / 4 // 45度
	g_est := 9.8

	// 从 y = x*tanθ - (g*x²)/(2*v₀²*cos²θ) 估算 v₀
	// 使用中间点估算
	x_mid := p2.X
	y_mid := p2.Y
	cos2 := math.Pow(math.Cos(theta_est), 2)
	v0_est := math.Sqrt((g_est * x_mid * x_mid) / (2 * (x_mid*math.Tan(theta_est) - y_mid) * cos2))

	// 确保正值
	if v0_est <= 0 || math.IsNaN(v0_est) {
		v0_est = 10.0 // 默认值
	}
	if g_est <= 0 {
		g_est = 9.8
	}

	return PhysicsParams{
		V0:    v0_est,
		Theta: theta_est,
		G:     g_est,
	}
}
