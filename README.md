# fft-zig

High-performance FFT (Fast Fourier Transform) implementation in Zig with SIMD optimizations.

## 当前状态 / Current Status

- ✅ **CI 测试通过** - All tests passing on Zig 0.14.1 and 0.15.1
- ✅ **FFT 精度正确** - FFT accuracy >99% for all test cases
- ✅ **支持大数据** - Supports large data FFT (up to 5M samples)

## 快速开始 / Quick Start

```bash
# 构建项目 / Build
zig build

# 运行测试 / Run tests
zig build test

# 格式化代码 / Format code
zig fmt build.zig src/
```

## 功能特性 / Features

- ✅ Radix-2 FFT (SIMD优化)
- ✅ Radix-4 FFT (SIMD优化)
- ✅ Mixed-radix FFT (Bluestein算法)
- ✅ 并行处理大数据 / Parallel processing for large datasets
- ✅ 实数到复数 FFT / Real-to-complex FFT
- ✅ 逆FFT / Inverse FFT

## 兼容性 / Compatibility

- Zig 0.14.1 ✅
- Zig 0.15.1 ✅

添加条件编译支持两个 Zig 版本的 ArrayList API 差异：

```zig
// Zig 0.14.1: allocator 存储在结构体中
var list = std.ArrayList(T).init(allocator);

// Zig 0.15.1: allocator 传递给方法
var list = std.ArrayList(T){ .items = &[_]T{}, .capacity = 0 };
list.append(allocator, item);
```

## 依赖关系图 / Dependency Graph

```
types.zig
   ↑
   │
 ┌───────────────┬───────────────┬───────────────┬───────────────┐
 │               │               │               │               │
fft_radix2.zig  fft_radix4.zig  utils.zig      base.zig
   ↑               ↑               ↑               ↑
   │               │               │               │
   └───────┬───────┴───────┬───────┴───────┬───────┘
           │               │               │
      fft_parallel.zig     │           ifft.zig
           ↑               │               ↑
           │               │               │
      fft_mixed.zig────────┘               │
           ↑                               │
           └───────────────┬───────────────┘
                           │
                      fft_r2c.zig
```

fft.zig → (types, utils, fft_radix2, fft_radix4, fft_parallel, fft_mixed, ifft, fft_r2c, base)

## 技术说明

### FFT 归一化约定

不同的 FFT 库有不同的归一化约定：

1. **无归一化**（本修复采用）
   - Forward FFT: X[k] = Σ x[n]·e^(-2πikn/N)
   - Inverse FFT: x[n] = (1/N)·Σ X[k]·e^(2πikn/N)
   - 优点：与数学定义一致，IFFT 恢复原值

2. **对称归一化**
   - Forward FFT: X[k] = (1/√N)·Σ x[n]·e^(-2πikn/N)
   - Inverse FFT: x[n] = (1/√N)·Σ X[k]·e^(2πikn/N)
   - 优点：对称，保持能量

3. **前向归一化**（原错误实现）
   - Forward FFT: X[k] = (1/N)·Σ x[n]·e^(-2πikn/N)
   - Inverse FFT: x[n] = Σ X[k]·e^(2πikn/N)
   - 缺点：不常见，IFFT 不能恢复原值

本实现采用第 1 种约定（无归一化），这是最常见的约定，被 FFTW、NumPy、MATLAB 等主流库采用。

### 验证方法

可以用以下方法验证 FFT 正确性：

1. **Parseval 定理**: Σ|x[n]|² = (1/N)·Σ|X[k]|²
2. **单位冲激**: FFT([1,0,0,...,0]) = [1,1,1,...,1]
3. **单频正弦**: FFT(sin(2πk₀n/N)) 在 k=k₀ 处有峰值 N/2

## 贡献 / Contributing

代码使用 `zig fmt` 格式化 / Code formatted with `zig fmt`

## 许可 / License

See repository license
