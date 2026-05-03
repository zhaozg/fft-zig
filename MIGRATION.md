# 迁移指南：从硬编码 f64 到泛型 FFT

本文档帮助你将现有代码从旧版（硬编码 `f64`）迁移到新版（泛型 `comptime T`）API。

## 破坏性变更

### 1. 函数签名变化

**旧版（v1.x）：**
```zig
pub fn fft(allocator: std.mem.Allocator, input: []const f64, output: []Complex) !void
pub fn fftInPlace(allocator: std.mem.Allocator, data: []Complex) !void
```

**新版（v2.x）：**
```zig
pub fn fft(comptime T: type, allocator: std.mem.Allocator, input: []const T, output: []std.math.Complex(T)) !void
pub fn fftInPlace(comptime T: type, allocator: std.mem.Allocator, data: []std.math.Complex(T)) !void
```

主要变化：
- 新增第一个参数 `comptime T: type`，用于指定浮点类型
- `Complex` 变为 `std.math.Complex(T)`
- 输入/输出类型由 `T` 决定

### 2. 类型变化

| 旧类型 | 新类型 |
|--------|--------|
| `Complex` | `std.math.Complex(T)` 或 `fft.Complex(T)` |
| `VectorF64` | `fft.Vector4(f64)` 或 `@Vector(4, T)` |
| `VectorF64x8` | `fft.Vector8(f64)` 或 `@Vector(8, T)` |
| `VectorComplex` | `fft.VectorComplex(f64)` 或自定义结构体 |

## 自动迁移步骤

### 步骤 1：添加泛型参数

在所有调用 `fft`、`fftInPlace`、`fftR2C` 的地方，在第一个参数位置添加类型参数：

```zig
// 旧版
try fft(allocator, input, output);

// 新版 - f64
try fft(f64, allocator, input, output);

// 新版 - f32
try fft(f32, allocator, input, output);
```

### 步骤 2：更新类型引用

```zig
// 旧版
const Complex = @import("fft").Complex;
var data: []Complex = ...;

// 新版
const fft = @import("fft");
var data: []std.math.Complex(f64) = ...;
// 或使用泛型辅助
var data_f32: []std.math.Complex(f32) = ...;
```

### 步骤 3：更新输入数据类型

```zig
// 旧版 - 只能使用 f64
const input: []const f64 = ...;

// 新版 - 可使用 f32 或 f64
const input_f32: []const f32 = ...;
const input_f64: []const f64 = ...;
```

## 向后兼容方案

如果你暂时不想修改现有代码，可以使用向后兼容别名：

```zig
const fft = @import("fft");

// 这些别名与旧版 API 签名相同
try fft.fft_f64(allocator, input, output);       // 等价于 fft(f64, ...)
try fft.fftInPlace_f64(allocator, data);          // 等价于 fftInPlace(f64, ...)
try fft.fftR2C_f64(allocator, input, output, mag); // 等价于 fftR2C(f64, ...)
```

## 新特性：使用 f32

迁移到泛型 API 后，你可以轻松使用 `f32` 获得更好的性能：

```zig
const std = @import("std");
const fft = @import("fft");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const N = 1024;

    // f32 输入
    const input = try allocator.alloc(f32, N);
    defer allocator.free(input);

    // ... 填充输入数据 ...

    const output = try allocator.alloc(std.math.Complex(f32), N);
    defer allocator.free(output);

    // 使用 f32 进行 FFT
    try fft.fft(f32, allocator, input, output);
}
```

## 性能对比

| 特性 | f32 | f64 |
|------|-----|-----|
| 内存占用 | 每个复数 8 字节 | 每个复数 16 字节 |
| 相对速度 | 1.5x - 2.5x 更快 | 基准 |
| 有效精度 | ~7 位小数 | ~15 位小数 |
| 适用场景 | 音频、图像、实时处理 | 科学计算、高精度分析 |

## 常见问题

### Q: 我可以在同一个程序中同时使用 f32 和 f64 吗？

**可以。** 泛型设计允许在同一程序中混合使用：

```zig
// f32 FFT
try fft.fft(f32, allocator, input_f32, output_f32);

// f64 FFT
try fft.fft(f64, allocator, input_f64, output_f64);
```

### Q: 旧版代码中的 `Complex` 类型还能用吗？

`fft.Complex` 仍然存在，但它是 `Complex(f64)` 的别名。建议迁移到 `std.math.Complex(T)`。

### Q: 迁移后性能会下降吗？

不会。泛型是编译期特性，`fft(f64, ...)` 生成的代码与旧版硬编码 `f64` 版本性能相同（±5% 以内）。

### Q: 如何验证迁移正确性？

运行测试套件：
```bash
zig build test
```

所有测试会同时验证 `f32` 和 `f64` 的正确性。
