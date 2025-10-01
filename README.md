# fft-zig

High-performance FFT (Fast Fourier Transform) implementation in Zig with SIMD optimizations.

## 当前状态 / Current Status

✅ **CI 测试通过** - All tests passing on Zig 0.14.1 and 0.15.1  
✅ **FFT 精度正确** - FFT accuracy >90% for all test cases  
✅ **支持大数据** - Supports large data FFT (up to 5M samples)  
⚠️ **Radix-4 暂时禁用** - Radix-4 FFT temporarily disabled (using radix-2 fallback)

### 已知问题 / Known Issues

Radix-4 FFT 实现存在 bug，导致输出幅值约为预期的 30.7%。当前使用 radix-2 作为替代方案，性能影响约 5-10%。

The radix-4 FFT implementation has a bug causing output magnitudes to be ~30.7% of expected. Currently using radix-2 as fallback with ~5-10% performance impact.

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
- ⚠️ Radix-4 FFT (暂时禁用/temporarily disabled)
- ✅ Mixed-radix FFT (Bluestein算法)
- ✅ 并行处理大数据 / Parallel processing for large datasets
- ✅ 实数到复数 FFT / Real-to-complex FFT
- ✅ 逆FFT / Inverse FFT

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

## 测试结果 / Test Results

- ✅ 25/25 tests passing on Zig 0.14.1
- ✅ 25/25 tests passing on Zig 0.15.1
- ✅ Large data validation (1M, 2M, 4M, 5M samples)
- ✅ FFT accuracy >90% for all test cases

## 兼容性 / Compatibility

- Zig 0.14.1 ✅
- Zig 0.15.1 ✅

支持两个版本的 ArrayList API 差异 / Supports ArrayList API differences between versions

## 贡献 / Contributing

代码使用 `zig fmt` 格式化 / Code formatted with `zig fmt`

## 许可 / License

See repository license
