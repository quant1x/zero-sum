# Fitting 包

此包提供非线性优化算法，使用Levenberg-Marquardt方法从观测点反推抛物线物理参数。

## 功能特性

- **SVD拟合**：使用奇异值分解进行初始抛物线拟合
- **LM优化**：迭代优化参数，包含阻尼和约束
- **QR求解**：高效解线性方程组

## 使用方法

```go
import "github.com/quant1x/zero-sum/fitting"

points := [3]fitting.Point{
    {X: 1, Y: 6.18},
    {X: 2, Y: 6.32},
    {X: 3, Y: 6.09},
}

solver := fitting.ProjectileSolver{
    Points:  points,
    Tol:     1e-6,
    MaxIter: 1000,
}

params, err := solver.Solve()
if err != nil {
    log.Fatal(err)
}

fmt.Printf("初速度: %.2f m/s\n", params.V0)
fmt.Printf("抛射角: %.2f°\n", params.Theta*180/math.Pi)
fmt.Printf("重力: %.2f m/s²\n", params.G)
```

## 算法流程

1. **SVD拟合**：三点拟合抛物线，得到初始a, b, c
2. **参数估计**：从抛物线参数推导物理参数
3. **LM优化**：迭代最小化残差，更新参数

## 测试

运行测试：

```bash
go test ./fitting
go test -bench=. ./fitting
```
