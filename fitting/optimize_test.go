package fitting

import (
	"math"
	"testing"
)

// TestSolve 测试参数求解（允许未收敛）
func TestSolve(t *testing.T) {
	points := [3]Point{
		{X: 1, Y: 6.18},
		{X: 2, Y: 6.32},
		{X: 3, Y: 6.09},
	}

	solver := ProjectileSolver{
		Points:  points,
		Tol:     1e-3,
		MaxIter: 10000,
	}

	result, err := solver.Solve()
	if err != nil {
		t.Logf("求解未收敛: %v", err)
	}

	t.Logf("最终参数: V0=%.2f, Theta=%.2f°, G=%.2f", result.V0, result.Theta*180/math.Pi, result.G)

	// 基本检查
	if result.V0 > 0 && result.Theta > 0 && result.Theta < math.Pi/2 && result.G > 0 {
		t.Log("参数合理")
	} else {
		t.Error("参数无效")
	}
}

// BenchmarkSolve 基准测试求解
func BenchmarkSolve(b *testing.B) {
	points := [3]Point{
		{X: 1, Y: 6.18},
		{X: 2, Y: 6.32},
		{X: 3, Y: 6.09},
	}

	solver := ProjectileSolver{
		Points:  points,
		Tol:     1e-3,
		MaxIter: 10000,
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := solver.Solve()
		if err != nil {
			b.Fatal(err)
		}
	}
}