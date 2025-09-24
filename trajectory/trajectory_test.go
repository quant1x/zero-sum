package trajectory

import (
	"fmt"
	"math"
	"testing"
)

// TestOptimizeLM 测试LM优化算法
func TestOptimizeLM(t *testing.T) {
	points := []Point{
		{X: 1, Y: 6.18},
		{X: 2, Y: 6.32},
		{X: 3, Y: 6.09},
	}

	initialParams := []float64{
		math.Atan(1.5), // 初始角度
		5.0,            // 初始速度
		9.8,            // 初始重力
	}

	params, err := OptimizeLM(points, initialParams, 100, 1e-6)
	if err != nil {
		t.Fatal(err)
	}

	theta, v0, g := params[0], params[1], params[2]
	fmt.Printf("优化结果: 角度=%.2f°, 初速度=%.2f m/s, 重力=%.2f m/s²\n", theta*180/math.Pi, v0, g)

	if theta <= 0 || theta >= math.Pi/2 {
		t.Errorf("角度超出范围: %f", theta)
	}
	if v0 <= 0 {
		t.Errorf("速度无效: %f", v0)
	}
	if g <= 0 {
		t.Errorf("重力无效: %f", g)
	}
}

// BenchmarkOptimizeLM 基准测试LM优化
func BenchmarkOptimizeLM(b *testing.B) {
	points := []Point{
		{X: 1, Y: 6.18},
		{X: 2, Y: 6.32},
		{X: 3, Y: 6.09},
	}

	initialParams := []float64{math.Atan(1.5), 5.0, 9.8}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := OptimizeLM(points, initialParams, 100, 1e-6)
		if err != nil {
			b.Fatal(err)
		}
	}
}