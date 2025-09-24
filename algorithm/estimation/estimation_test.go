package estimation

import (
	"math"
	"testing"
)

func TestProjectileSolver_Solve(t *testing.T) {
	// 使用能形成抛物线的合理测试数据
	// 基于参数：v0=10, θ=45°, g=9.8
	trueParams := PhysicsParams{V0: 10, Theta: math.Pi / 4, G: 9.8}

	points := [3]Point{}
	for i, x := range []float64{2, 5, 8} {
		y := x*math.Tan(trueParams.Theta) - (trueParams.G*x*x)/(2*trueParams.V0*trueParams.V0*math.Pow(math.Cos(trueParams.Theta), 2))
		points[i] = Point{X: x, Y: y}
	}

	solver := ProjectileSolver{
		Points:  points,
		Tol:     1e-6,
		MaxIter: 1000, // 增加迭代次数
	}

	initial := EstimateInitial(points)
	result, err := solver.Solve(initial)
	if err != nil {
		t.Fatalf("求解失败: %v", err)
	}

	// 检查结果是否合理（正值且在合理范围内）
	if result.V0 <= 0 || result.G <= 0 || result.Theta < 0 || result.Theta > math.Pi/2 {
		t.Errorf("参数超出合理范围: %v", result)
	}

	// 计算残差检查拟合质量
	residuals := solver.residuals(result)
	totalError := norm(residuals)
	if totalError > 0.1 { // 允许小误差
		t.Errorf("拟合误差过大: %f", totalError)
	}

	t.Logf("真实参数: V0=%.3f, Theta=%.3f, G=%.3f", trueParams.V0, trueParams.Theta, trueParams.G)
	t.Logf("估算参数: V0=%.3f, Theta=%.3f, G=%.3f", result.V0, result.Theta, result.G)
}

func TestEstimateInitial(t *testing.T) {
	points := [3]Point{
		{X: 0, Y: 0},
		{X: 5, Y: 2.5},
		{X: 10, Y: 0},
	}

	initial := EstimateInitial(points)
	if initial.V0 <= 0 || initial.G <= 0 {
		t.Error("初始估算参数无效")
	}
}

func BenchmarkProjectileSolver_Solve(b *testing.B) {
	points := [3]Point{
		{X: 2, Y: 1.414},
		{X: 4, Y: 2.828},
		{X: 6, Y: 4.242},
	}

	solver := ProjectileSolver{
		Points:  points,
		Tol:     1e-6,
		MaxIter: 100,
	}

	initial := PhysicsParams{V0: 10, Theta: math.Pi / 4, G: 9.8}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = solver.Solve(initial)
	}
}
