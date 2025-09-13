package main

import (
	"fmt"
	"math"
)

type Point struct{ X, Y float64 }

type ProjectileSolver struct {
	Points  [3]Point
	Tol     float64
	MaxIter int
}

type PhysicsParams struct {
	V0, Theta, G float64
}

// 改进的求解方法（Levenberg-Marquardt 算法）
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
			return p, fmt.Errorf("矩阵奇异")
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
	return p, fmt.Errorf("超过最大迭代次数")
}

// 计算残差向量
func (s *ProjectileSolver) residuals(p PhysicsParams) [3]float64 {
	var F [3]float64
	for i, pt := range s.Points {
		x, y := pt.X, pt.Y
		F[i] = y - (x*math.Tan(p.Theta) - (p.G*x*x)/(2*p.V0*p.V0*math.Pow(math.Cos(p.Theta), 2)))
	}
	return F
}

// 计算雅可比矩阵（数值微分法）
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

// 矩阵运算工具函数
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

func v1solveLinearSystem(A [3][3]float64, b [3]float64) ([3]float64, bool) {
	// 此处使用高斯消元法（实际应用建议使用优化库）
	// 此处简化为调用前文实现的 solveLinearSystem 函数
	// 可替换为更稳健的线性代数库实现
	return [3]float64{}, true
}

func norm(v [3]float64) float64 {
	return math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}

// 自动估算初始值
func estimateInitial(points [3]Point) PhysicsParams {
	// 基于抛物线顶点估算
	x1, y1 := points[0].X, points[0].Y
	x2, y2 := points[1].X, points[1].Y
	x3, y3 := points[2].X, points[2].Y

	// 假设对称轴位于中间点
	xVertex := (x1 + x3) / 2
	gEstimate := 9.8 // 初始假设
	v0Estimate := math.Sqrt((x3 * x3 * gEstimate) / (2 * (x3*math.Tan(math.Pi/4) - y3)))

	_ = y1
	_ = x2
	_ = y2
	_ = xVertex
	return PhysicsParams{
		V0:    v0Estimate,
		Theta: math.Pi / 4,
		G:     gEstimate,
	}
}

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

func main() {
	// 生成测试数据（已知参数：v0=20m/s, θ=45°, g=9.81m/s²）
	trueParams := PhysicsParams{
		V0:    20,
		Theta: math.Pi / 4,
		G:     9.81,
	}

	// 计算三个理论轨迹点
	//points := [3]Point{
	//	{X: 5, Y: 5*1 - 9.81*5*5/(2*20*20*0.5)},
	//	{X: 10, Y: 10*1 - 9.81*10*10/(2*20*20*0.5)},
	//	{X: 15, Y: 15*1 - 9.81*15*15/(2*20*20*0.5)},
	//}
	// 更合理的测试数据（三点不共线）
	points := [3]Point{
		{X: 2.0, Y: 5.0},
		{X: 5.0, Y: 8.0},
		{X: 8.0, Y: 5.0},
	}

	solver := ProjectileSolver{
		Points:  points,
		Tol:     1e-6,
		MaxIter: 100,
	}

	// 初始猜测（带误差的猜测值）
	//initial := PhysicsParams{
	//	V0:    18,
	//	Theta: math.Pi / 3.5,
	//	G:     8.5,
	//}

	// 自动估算初始值
	initial := estimateInitial(points)

	// 执行求解
	result, err := solver.Solve(initial)
	if err != nil {
		fmt.Println("求解失败:", err)
		return
	}

	// 输出结果
	fmt.Println("=== 真实参数 ===")
	fmt.Printf("初速度: %.2f m/s\n", trueParams.V0)
	fmt.Printf("抛射角: %.1f°\n", trueParams.Theta*180/math.Pi)
	fmt.Printf("重力加速度: %.2f m/s²\n\n", trueParams.G)

	fmt.Println("=== 反推结果 ===")
	fmt.Printf("初速度: %.2f m/s (误差: %.2f%%)\n",
		result.V0, 100*math.Abs((result.V0-trueParams.V0)/trueParams.V0))
	fmt.Printf("抛射角: %.1f° (误差: %.2f%%)\n",
		result.Theta*180/math.Pi, 100*math.Abs((result.Theta-trueParams.Theta)/trueParams.Theta))
	fmt.Printf("重力加速度: %.2f m/s² (误差: %.2f%%)",
		result.G, 100*math.Abs((result.G-trueParams.G)/trueParams.G))
}
