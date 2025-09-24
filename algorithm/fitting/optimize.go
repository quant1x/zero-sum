// Package fitting 提供非线性优化算法，使用Levenberg-Marquardt方法反推抛物线物理参数。
package fitting

import (
	"errors"
	"fmt"
	"gonum.org/v1/gonum/mat"
	"log"
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

// Solve 求解抛物线参数。
func (s *ProjectileSolver) Solve() (PhysicsParams, error) {
	a, b, _, ok := s.fitParabola()
	if !ok {
		return PhysicsParams{}, errors.New("三点不构成抛物线")
	}

	initial := s.estimateFromParabola(a, b)
	return s.refineWithLM(initial)
}

// fitParabola 三点抛物线拟合 y = ax² + bx + c
func (s *ProjectileSolver) fitParabola() (a, b, c float64, ok bool) {
	x1, y1 := s.Points[0].X, s.Points[0].Y
	x2, y2 := s.Points[1].X, s.Points[1].Y
	x3, y3 := s.Points[2].X, s.Points[2].Y

	A := mat.NewDense(3, 3, []float64{
		x1 * x1, x1, 1,
		x2 * x2, x2, 1,
		x3 * x3, x3, 1,
	})
	B := mat.NewVecDense(3, []float64{y1, y2, y3})

	var svd mat.SVD
	if !svd.Factorize(A, mat.SVDThin) {
		return 0, 0, 0, false
	}

	const rcond = 1e-15
	rank := svd.Rank(rcond)
	if rank == 0 {
		log.Fatal("zero rank system")
	}
	var x mat.VecDense
	svd.SolveVecTo(&x, B, rank)

	return x.AtVec(0), x.AtVec(1), x.AtVec(2), true
}

// estimateFromParabola 从抛物线参数估计物理参数
func (s *ProjectileSolver) estimateFromParabola(a, b float64) PhysicsParams {
	theta := math.Atan(b)
	cos := math.Cos(theta)

	// 使用射程公式估算初速度
	maxX := s.findMaxX()
	v0 := math.Sqrt((maxX * 9.8) / math.Sin(2*theta))

	// 通过抛物线系数计算重力加速度
	g := math.Abs(-2 * a * v0 * v0 * cos * cos)

	return PhysicsParams{
		V0:    math.Max(v0, 1.0),
		Theta: clamp(theta, 1e-3, math.Pi/2-1e-3),
		G:     math.Max(g, 1.0),
	}
}

// refineWithLM 带阻尼的最小二乘优化（Levenberg-Marquardt）
func (s *ProjectileSolver) refineWithLM(initial PhysicsParams) (PhysicsParams, error) {
	p := initial
	lambda := 0.1
	adjustFactor := 2.0
	prevErr := math.MaxFloat64

	for iter := 0; iter < s.MaxIter; iter++ {
		res := s.residuals(p)
		J := s.numericalJacobian(p)

		// 构建线性系统 (J^T J + λI)Δ = J^T r
		var JTJ mat.Dense
		JTJ.Mul(J.T(), J)
		addDiagonal(&JTJ, lambda)

		JTr := mat.NewVecDense(3, nil)
		rVec := mat.NewVecDense(len(res), res)
		JTr.MulVec(J.T(), rVec)

		// 解线性系统
		delta, ok := solveLinearSystem(&JTJ, JTr)
		if !ok {
			lambda *= adjustFactor
			continue
		}

		// 试验参数更新
		trial := PhysicsParams{
			V0:    p.V0 - delta.At(0, 0),
			Theta: p.Theta - delta.At(1, 0),
			G:     p.G - delta.At(2, 0),
		}

		// 物理约束
		trial.V0 = math.Max(trial.V0, 0.1)
		trial.G = math.Max(trial.G, 0.1)
		trial.Theta = clamp(trial.Theta, 1e-5, math.Pi/2-1e-5)

		// 计算新残差
		newRes := s.residuals(trial)
		newNorm := mat.Norm(mat.NewVecDense(len(newRes), newRes), 2)

		if newNorm < prevErr { // 接受更新
			p = trial
			prevErr = newNorm
			lambda /= adjustFactor
			adjustFactor = math.Max(adjustFactor*0.9, 1.1)

			if newNorm < s.Tol {
				return p, nil
			}
		} else { // 拒绝更新
			lambda *= adjustFactor
			adjustFactor = math.Min(adjustFactor*1.1, 10.0)
		}
	}
	return p, fmt.Errorf("超过最大迭代次数")
}

// numericalJacobian 数值方法计算雅可比矩阵
func (s *ProjectileSolver) numericalJacobian(p PhysicsParams) *mat.Dense {
	epsilon := 1e-8
	res0 := s.residuals(p)
	J := mat.NewDense(3, 3, nil)

	for col := 0; col < 3; col++ {
		// 创建参数扰动
		perturbed := p
		switch col {
		case 0:
			perturbed.V0 += epsilon
		case 1:
			perturbed.Theta += epsilon
		case 2:
			perturbed.G += epsilon
		}

		resPerturbed := s.residuals(perturbed)
		for row := 0; row < 3; row++ {
			deriv := (resPerturbed[row] - res0[row]) / epsilon
			J.Set(row, col, deriv)
		}
	}
	return J
}

// residuals 计算残差向量
func (s *ProjectileSolver) residuals(p PhysicsParams) []float64 {
	res := make([]float64, 3)
	for i, pt := range s.Points {
		term1 := pt.X * math.Tan(p.Theta)
		term2 := (p.G * pt.X * pt.X) / (2 * p.V0 * p.V0 * math.Pow(math.Cos(p.Theta), 2))
		res[i] = pt.Y - (term1 - term2)
	}
	return res
}

// addDiagonal 矩阵对角线增强
func addDiagonal(m *mat.Dense, lambda float64) {
	r, c := m.Dims()
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			if i == j {
				m.Set(i, j, m.At(i, j)*(1+lambda))
			}
		}
	}
}

// solveLinearSystem 求解线性方程组（使用QR分解）
func solveLinearSystem(A *mat.Dense, b *mat.VecDense) (*mat.Dense, bool) {
	var qr mat.QR
	qr.Factorize(A)

	var x mat.Dense
	if err := x.Solve(&qr, b); err != nil {
		return nil, false
	}
	return &x, true
}

// findMaxX 辅助函数
func (s *ProjectileSolver) findMaxX() float64 {
	max := 0.0
	for _, p := range s.Points {
		if p.X > max {
			max = p.X
		}
	}
	return max
}

func clamp(value, min, max float64) float64 {
	return math.Max(min, math.Min(value, max))
}