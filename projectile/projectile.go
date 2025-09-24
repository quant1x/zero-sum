// Package projectile 提供斜抛运动模拟，计算轨迹参数和生成运动路径。
package projectile

import (
	"math"
)

// Projectile 表示斜抛运动参数。
type Projectile struct {
	V0      float64 // 初始速度 (m/s)
	Angle   float64 // 抛射角度 (度)
	Gravity float64 // 重力加速度 (m/s²)
}

// TrajectoryPoint 表示轨迹上的一个点。
type TrajectoryPoint struct {
	Time float64 // 时间 (秒)
	X    float64 // 水平位置 (米)
	Y    float64 // 垂直位置 (米)
	Vx   float64 // 水平速度 (m/s)
	Vy   float64 // 垂直速度 (m/s)
}

// Calculate 计算斜抛运动的关键参数。
// 返回：最大高度、总飞行时间、最大射程。
func (p *Projectile) Calculate() (maxHeight, totalTime, maxRange float64) {
	theta := p.Angle * math.Pi / 180 // 角度转弧度
	v0x := p.V0 * math.Cos(theta)
	v0y := p.V0 * math.Sin(theta)

	// 最高点时间
	tPeak := v0y / p.Gravity
	// 总飞行时间
	totalTime = 2 * tPeak
	// 最大高度
	maxHeight = v0y*tPeak - 0.5*p.Gravity*tPeak*tPeak
	// 最大射程
	maxRange = v0x * totalTime

	return maxHeight, totalTime, maxRange
}

// GenerateTrajectory 生成轨迹数据。
// steps: 时间步数，返回完整的轨迹点数组。
func (p *Projectile) GenerateTrajectory(steps int) []TrajectoryPoint {
	theta := p.Angle * math.Pi / 180
	v0x := p.V0 * math.Cos(theta)
	v0y := p.V0 * math.Sin(theta)

	_, totalTime, _ := p.Calculate()
	dt := totalTime / float64(steps)

	var trajectory []TrajectoryPoint
	for t := 0.0; t <= totalTime; t += dt {
		x := v0x * t
		y := v0y*t - 0.5*p.Gravity*t*t
		vy := v0y - p.Gravity*t
		trajectory = append(trajectory, TrajectoryPoint{
			Time: t,
			X:    x,
			Y:    y,
			Vx:   v0x,
			Vy:   vy,
		})
	}
	return trajectory
}