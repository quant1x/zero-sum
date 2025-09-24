package regression

import (
	"testing"
)

// BenchmarkQuadraticRegression 测试二次回归的性能
func BenchmarkQuadraticRegression(b *testing.B) {
	// 准备测试数据：10个点
	points := []Point{
		{0, 5.64}, {1, 6.18}, {2, 6.32}, {3, 6.09},
		{4, 5.47}, {5, 4.50}, {6, 3.18}, {7, 1.51},
		{8, -0.61}, {9, -3.18},
	}

	b.ResetTimer() // 重置计时器，排除初始化时间
	for i := 0; i < b.N; i++ {
		_, _, err := QuadraticRegression(points)
		if err != nil {
			b.Fatal(err)
		}
	}
}

// BenchmarkQuadraticRegressionLarge 测试更大数据集的性能
func BenchmarkQuadraticRegressionLarge(b *testing.B) {
	// 准备更大测试数据：100个点
	points := make([]Point, 100)
	for i := 0; i < 100; i++ {
		x := float64(i)
		y := 0.5*x*x - 2*x + 10 + float64(i%10) // 添加一些噪声
		points[i] = Point{X: x, Y: y}
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _, err := QuadraticRegression(points)
		if err != nil {
			b.Fatal(err)
		}
	}
}

// BenchmarkSolveGauss 测试高斯消元法的性能
func BenchmarkSolveGauss(b *testing.B) {
	// 3x3矩阵（对应二次回归）
	A := [][]float64{
		{1, 1, 1},
		{4, 2, 1},
		{9, 3, 1},
	}
	y := []float64{1, 2, 3}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := solveGauss(A, y)
		if err != nil {
			b.Fatal(err)
		}
	}
}