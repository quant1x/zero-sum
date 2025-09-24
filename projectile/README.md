# Projectile 包

此包提供斜抛运动模拟，计算关键参数并生成完整运动轨迹。

## 功能特性

- **参数计算**：最大高度、飞行时间、射程
- **轨迹生成**：离散化时间点，包含位置和速度
- **物理准确**：基于经典力学公式

## 使用方法

```go
import "github.com/quant1x/zero-sum/projectile"

p := projectile.Projectile{
    V0:      20,   // 初速度 m/s
    Angle:   45,   // 角度 度
    Gravity: 9.81, // 重力 m/s²
}

// 计算参数
maxHeight, totalTime, maxRange := p.Calculate()

// 生成轨迹
trajectory := p.GenerateTrajectory(100) // 100个时间点
for _, point := range trajectory {
    fmt.Printf("t=%.2f, x=%.2f, y=%.2f\n", point.Time, point.X, point.Y)
}
```

## 公式

- 位置：`x = v₀·cos(θ)·t`, `y = v₀·sin(θ)·t - 0.5·g·t²`
- 速度：`vx = v₀·cos(θ)`, `vy = v₀·sin(θ) - g·t`

## 测试

运行测试：

```bash
go test ./projectile
go test -bench=. ./projectile
```
