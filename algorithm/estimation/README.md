# estimation

参数估计算法包，使用Levenberg-Marquardt方法反推抛物线物理参数。

## 功能特性

- **Levenberg-Marquardt算法**：非线性最小二乘优化，用于参数估计
- **数值雅可比矩阵**：自动计算偏导数，无需手动推导
- **自动初始值估算**：基于观测点智能估算初始参数
- **物理约束**：确保参数在合理范围内（正值、角度限制）
- **高性能**：纯Go实现，无外部依赖

## 使用方法

```go
package main

import (
    "fmt"
    "log"
    "math"

    "github.com/your-repo/estimation"
)

func main() {
    // 定义观测点
    points := [3]estimation.Point{
        {X: 2.0, Y: 1.414},
        {X: 4.0, Y: 2.828},
        {X: 6.0, Y: 4.242},
    }

    // 创建求解器
    solver := estimation.ProjectileSolver{
        Points:  points,
        Tol:     1e-6,    // 收敛容差
        MaxIter: 100,     // 最大迭代次数
    }

    // 估算初始值
    initial := estimation.EstimateInitial(points)

    // 执行优化
    result, err := solver.Solve(initial)
    if err != nil {
        log.Fatal(err)
    }

    fmt.Printf("估算参数: V0=%.3f, Theta=%.3f rad, G=%.3f\n",
        result.V0, result.Theta, result.G)
}
```

## API文档

### 类型定义

#### Point
```go
type Point struct {
    X, Y float64
}
```
二维坐标点。

#### PhysicsParams
```go
type PhysicsParams struct {
    V0    float64 // 初速度 (m/s)
    Theta float64 // 抛射角 (弧度)
    G     float64 // 重力加速度 (m/s²)
}
```
抛物线运动物理参数。

#### ProjectileSolver
```go
type ProjectileSolver struct {
    Points  [3]Point // 观测点
    Tol     float64  // 收敛容差
    MaxIter int      // 最大迭代次数
}
```
参数求解器。

### 函数

#### Solve
```go
func (s *ProjectileSolver) Solve(initial PhysicsParams) (PhysicsParams, error)
```
执行Levenberg-Marquardt优化，返回估算参数。

#### EstimateInitial
```go
func EstimateInitial(points [3]Point) PhysicsParams
```
基于观测点自动估算初始参数值。

## 算法原理

该包使用Levenberg-Marquardt算法最小化观测点与理论抛物线轨迹之间的残差：

```
y = x·tan(θ) - (g·x²)/(2·v₀²·cos²(θ))
```

其中：
- `v₀`: 初速度
- `θ`: 抛射角
- `g`: 重力加速度

算法特点：
- **数值微分**：使用有限差分计算雅可比矩阵
- **阻尼调节**：动态调整阻尼系数以提高收敛性
- **线性求解**：使用伴随矩阵法解正规方程

## 性能基准

```
BenchmarkProjectileSolver_Solve-8   	  100000	     10345 ns/op
```

## 依赖关系

无外部依赖，纯标准库实现。

## 许可证

MIT License