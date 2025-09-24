# Regression 包

此包提供二次回归分析，用于将数据点拟合到二次多项式 y = ax² + bx + c 并分析曲线属性。

## 功能特性

- **二次回归**：将点拟合到 y = ax² + bx + c
- **预测**：为新 x 值生成预测函数
- **运动分析**：分析曲线开口和运动方向

## 使用方法

```go
import "github.com/quant1x/zero-sum/regression"

points := []regression.Point{
    {0, 5.64}, {1, 6.18}, {2, 6.32}, {3, 6.09},
}

predict, coeffs, err := regression.QuadraticRegression(points)
if err != nil {
    log.Fatal(err)
}

fmt.Printf("方程: y = %.2fx² + %.2fx + %.2f\n", coeffs[0], coeffs[1], coeffs[2])

// 预测新值
y := predict(4.0)

// 分析运动
direction, opening := regression.AnalyzeMotion(coeffs[0], coeffs[1], coeffs[2], 4.0)
```

## 要求

- 至少 3 个数据点
- 点不应共线（避免奇异矩阵）

## 性能

- 小数据集（<100 个点）：快速，适合实时分析
- 大数据集：矩阵运算 O(n)，生产环境考虑优化

## 测试

运行基准测试：

```bash
go test -bench=. ./regression
```

