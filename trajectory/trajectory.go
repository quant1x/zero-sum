// Package trajectory 提供抛物线轨迹参数优化，使用Levenberg-Marquardt算法反推物理参数。
package trajectory

import (
	"errors"
	"math"
)

// Point 表示具有 X 和 Y 坐标的二维点。
type Point struct{ X, Y float64 }

// Residuals 计算目标函数残差（抛物线模型）。
func Residuals(params []float64, points []Point, res []float64) {
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

// Jacobian 计算雅可比矩阵。
func Jacobian(params []float64, points []Point, jac [][]float64) {
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

// OptimizeLM 使用LM算法优化抛物线参数。
// 输入：点集、初始参数、最大迭代次数、收敛容差。
// 输出：优化后的参数 [θ, v0, g] 和错误。
func OptimizeLM(points []Point, initialParams []float64, maxIter int, tol float64) ([]float64, error) {
	if len(points) < 3 {
		return nil, errors.New("需要至少3个点")
	}
	params := make([]float64, len(initialParams))
	copy(params, initialParams)

	lambda := 0.01
	res := make([]float64, len(points))
	jac := make([][]float64, len(points))
	for i := range jac {
		jac[i] = make([]float64, 3)
	}

	for iter := 0; iter < maxIter; iter++ {
		Residuals(params, points, res)
		Jacobian(params, points, jac)

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
			JTJ[i][i] += lambda
		}

		// 求解增量方程
		delta := make([]float64, 3)
		if err := solve3x3(JTJ, JTr, delta); err != nil {
			return nil, err
		}

		// 尝试更新参数
		newParams := make([]float64, 3)
		copy(newParams, params)
		for i := 0; i < 3; i++ {
			newParams[i] -= delta[i] // 注意符号
		}

		// 物理约束
		newParams[0] = math.Max(1e-5, math.Min(newParams[0], math.Pi/2-1e-5)) // θ
		newParams[1] = math.Max(0.1, newParams[1])                           // v0
		newParams[2] = math.Max(0.1, newParams[2])                           // g

		// 计算残差变化
		newRes := make([]float64, len(points))
		Residuals(newParams, points, newRes)
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

// solve3x3 高斯消元法解3x3线性方程组。
func solve3x3(A [][]float64, b []float64, x []float64) error {
	// 增广矩阵
	aug := make([][]float64, 3)
	for i := 0; i < 3; i++ {
		aug[i] = make([]float64, 4)
		copy(aug[i][:3], A[i])
		aug[i][3] = b[i]
	}

	// 高斯消元
	for col := 0; col < 3; col++ {
		// 选主元
		maxRow := col
		for r := col; r < 3; r++ {
			if math.Abs(aug[r][col]) > math.Abs(aug[maxRow][col]) {
				maxRow = r
			}
		}
		aug[col], aug[maxRow] = aug[maxRow], aug[col]

		pivot := aug[col][col]
		if math.Abs(pivot) < 1e-12 {
			return errors.New("矩阵奇异")
		}

		// 归一化
		for j := col; j < 4; j++ {
			aug[col][j] /= pivot
		}

		// 消元
		for r := 0; r < 3; r++ {
			if r != col {
				factor := aug[r][col]
				for j := col; j < 4; j++ {
					aug[r][j] -= factor * aug[col][j]
				}
			}
		}
	}

	for i := 0; i < 3; i++ {
		x[i] = aug[i][3]
	}
	return nil
}

func norm(v []float64) float64 {
	sum := 0.0
	for _, val := range v {
		sum += val * val
	}
	return math.Sqrt(sum)
}