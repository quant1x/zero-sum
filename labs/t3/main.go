package main

import (
	"fmt"
	"math"
)

// 斜抛运动参数结构体
type Projectile struct {
	V0      float64 // 初始速度 (m/s)
	Angle   float64 // 抛射角度 (度)
	Gravity float64 // 重力加速度 (m/s²)
}

// 轨迹点结构体
type TrajectoryPoint struct {
	Time float64 // 时间 (秒)
	X    float64 // 水平位置 (米)
	Y    float64 // 垂直位置 (米)
	Vx   float64 // 水平速度 (m/s)
	Vy   float64 // 垂直速度 (m/s)
}

// 计算斜抛运动参数
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

// 生成轨迹数据
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

// 打印轨迹表格
func PrintTrajectory(trajectory []TrajectoryPoint) {
	fmt.Printf("%-8s %-8s %-8s %-8s %-8s\n", "Time", "X", "Y", "Vx", "Vy")
	for _, p := range trajectory {
		fmt.Printf("%-8.2f %-8.2f %-8.2f %-8.2f %-8.2f\n",
			p.Time, p.X, p.Y, p.Vx, p.Vy)
	}
}

func main() {
	// 测试案例：初速度 20m/s，角度 45度
	projectile := Projectile{
		V0:      20,
		Angle:   45,
		Gravity: 9.81,
	}

	projectile = Projectile{
		V0:      5.00,
		Angle:   56.31,
		Gravity: 9.80,
	}

	// 计算关键参数
	maxHeight, totalTime, maxRange := projectile.Calculate()
	fmt.Println("=== 斜抛运动关键参数 ===")
	fmt.Printf("最大高度: %.2f 米\n", maxHeight)
	fmt.Printf("总飞行时间: %.2f 秒\n", totalTime)
	fmt.Printf("水平射程: %.2f 米\n", maxRange)
	fmt.Println("\n=== 数值验证 ===")
	fmt.Printf("落地时 Y 坐标: %.2f 米\n",
		projectile.GenerateTrajectory(1)[1].Y) // 应接近 0

	// 生成轨迹数据（20个时间点）
	fmt.Println("\n=== 运动轨迹（前5个时间点）===")
	trajectory := projectile.GenerateTrajectory(20)
	PrintTrajectory(trajectory[:])
}
