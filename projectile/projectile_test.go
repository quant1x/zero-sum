package projectile

import (
	"math"
	"testing"
)

// TestCalculate 测试参数计算
func TestCalculate(t *testing.T) {
	p := Projectile{V0: 20, Angle: 45, Gravity: 9.81}

	maxHeight, totalTime, maxRange := p.Calculate()

	// 理论值计算
	theta := 45 * math.Pi / 180
	v0y := 20 * math.Sin(theta)
	tPeak := v0y / 9.81
	expectedHeight := v0y*tPeak - 0.5*9.81*tPeak*tPeak
	expectedTime := 2 * tPeak
	expectedRange := 20 * math.Cos(theta) * expectedTime

	if math.Abs(maxHeight-expectedHeight) > 1e-6 {
		t.Errorf("最大高度错误: 期望 %.6f, 得到 %.6f", expectedHeight, maxHeight)
	}
	if math.Abs(totalTime-expectedTime) > 1e-6 {
		t.Errorf("总时间错误: 期望 %.6f, 得到 %.6f", expectedTime, totalTime)
	}
	if math.Abs(maxRange-expectedRange) > 1e-6 {
		t.Errorf("射程错误: 期望 %.6f, 得到 %.6f", expectedRange, maxRange)
	}
}

// TestGenerateTrajectory 测试轨迹生成
func TestGenerateTrajectory(t *testing.T) {
	p := Projectile{V0: 5.00, Angle: 56.31, Gravity: 9.80}

	trajectory := p.GenerateTrajectory(10)

	if len(trajectory) != 11 { // 包括起始点
		t.Errorf("轨迹点数错误: 期望 11, 得到 %d", len(trajectory))
	}

	// 检查起始点
	if trajectory[0].Time != 0 || trajectory[0].X != 0 || trajectory[0].Y != 0 {
		t.Error("起始点错误")
	}

	// 检查落地点 (Y ≈ 0)
	last := trajectory[len(trajectory)-1]
	if math.Abs(last.Y) > 1e-3 {
		t.Errorf("落地点Y坐标错误: %.6f", last.Y)
	}
}

// BenchmarkGenerateTrajectory 基准测试轨迹生成
func BenchmarkGenerateTrajectory(b *testing.B) {
	p := Projectile{V0: 20, Angle: 45, Gravity: 9.81}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p.GenerateTrajectory(100)
	}
}