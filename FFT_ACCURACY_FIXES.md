# FFT 算法修复历史记录

> ⚠️ **注意**: 此文档记录了早期的 FFT 归一化修复。当前项目状态请参阅 **PLEASE_READ_FIRST.md**

## 历史背景

本文档记录了项目早期发现并修复的 FFT 归一化问题。这些修复已经完成并合并。

## 当前项目状态

**当前主要问题**: Radix-4 FFT 实现存在 bug（输出约为预期的 30.7%），已临时禁用，使用 radix-2 替代。

详情请参阅:
- **PLEASE_READ_FIRST.md** - 当前项目状态
- **README.md** - 项目概览
- **PR_SUMMARY.md** - 完整技术说明

---

## 历史修复记录 (已完成)

### 问题概述

原始问题：FFT 实现在数据处理时结果漂移严重，误差很大。

### 根本原因分析

#### 1. 前向 FFT 错误归一化（主要问题）

**位置**: `src/fft/fft_radix2.zig`

**问题**: 
- `fftRadix2()` 函数对前向 FFT 结果除以 N
- `fftRadix2SIMD()` 函数也有同样的问题
- 这导致 FFT 输出幅值缩小了 N 倍，造成严重的结果漂移

**原因**: 
标准 FFT 约定是：
- **前向 FFT**: 不归一化（输出幅值与 N 成正比）
- **逆向 IFFT**: 除以 N（恢复原始数值）

错误的归一化导致：
- FFT 结果比预期小 N 倍
- 如果再进行 IFFT，会再次除以 N，导致最终结果是原值的 1/N²
- 精度损失严重

**修复**: 已移除前向 FFT 中的归一化代码

#### 2. 缺失输出转换

**位置**: `src/fft/fft_r2c.zig`

**问题**: 对于中等大小输入，`fftR2C` 没有正确转换输出

**修复**: 已添加输出转换调用

#### 3. Bluestein 算法 chirp 镜像问题

**位置**: `src/fft/fft_mixed.zig`

**问题**: chirp 数组未正确镜像

**修复**: 已添加正确的镜像逻辑

#### 4. 测试验证不足

**位置**: `src/fft/fft_r2c.zig` 大数据测试

**问题**: 测试验证过于宽松

**修复**: 已加强测试验证标准

## 修复影响

### 精度改善
- FFT 输出现在正确缩放，与输入成比例
- 消除了 N 倍的幅值漂移
- FFT -> IFFT 循环现在能正确恢复原始数据

### 性能影响
- 移除归一化循环实际上略微提升了性能

### 测试覆盖
- 更新的测试现在能捕获精度回归
- 测试期望值与标准 FFT 行为一致

## 验证方法

可以用以下方法验证 FFT 正确性：

1. **Parseval 定理**: Σ|x[n]|² = (1/N)·Σ|X[k]|²
2. **单位冲激**: FFT([1,0,0,...,0]) = [1,1,1,...,1]
3. **单频正弦**: FFT(sin(2πk₀n/N)) 在 k=k₀ 处有峰值 N/2

修复后的实现满足所有这些性质。

---

*此文档保留用于历史参考。当前项目状态和待解决问题请参阅 PLEASE_READ_FIRST.md*
```zig
try fftInPlace(allocator, complex_buffer);
// 添加输出转换
convertToOutputSIMD(complex_buffer[0..out_len], output, magnitude);
```

### 3. Bluestein 算法 chirp 镜像问题

**位置**: `src/fft/fft_mixed.zig`

**问题**:
- Bluestein 算法用于非 2 的幂次长度 FFT
- chirp 数组 `b` 需要镜像以实现循环卷积
- 原实现只是零填充，导致非 2 的幂次 FFT 精度下降

**修复**:
```zig
// 修复前
for (n..m) |k| {
    b[k] = Complex{ .re = 0.0, .im = 0.0 };
}

// 修复后
// Mirror chirp values for circular convolution
for (1..n) |k| {
    b[m - k] = b[k];
}
for (n..(m - n + 1)) |k| {
    b[k] = Complex{ .re = 0.0, .im = 0.0 };
}
```

### 4. 测试验证不足

**位置**: `src/fft/fft_r2c.zig` 大数据测试

**问题**:
- 测试只验证 `magnitude[5] > 0.4`，过于宽松
- 对于单位幅值正弦波，频率 bin 的幅值应该约为 N/2

**修复**:
```zig
// 修复前
try expect(magnitude[5] > 0.4);

// 修复后
const expected_magnitude = @as(f64, @floatFromInt(size)) / 2.0;
try expect(magnitude[5] > expected_magnitude * 0.9); // 允许 10% 误差
```

## 修复影响

### 精度改善
- FFT 输出现在正确缩放，与输入成比例
- 消除了 N 倍的幅值漂移
- FFT -> IFFT 循环现在能正确恢复原始数据

### 性能影响
- 移除归一化循环实际上略微提升了性能
- 其他修复对性能影响可忽略

### 测试覆盖
- 更新的测试现在能捕获精度回归
- 测试期望值与标准 FFT 行为一致

## 修改文件列表

1. **src/fft/fft_radix2.zig**
   - 移除 `fftRadix2()` 的归一化
   - 移除 `fftRadix2SIMD()` 的归一化
   - 更新测试期望值

2. **src/fft/fft_r2c.zig**
   - 添加缺失的输出转换
   - 加强大数据测试验证

3. **src/fft/fft_mixed.zig**
   - 修复 Bluestein chirp 镜像
   - 更新测试参考 DFT 以匹配非归一化行为

## 验证步骤

需要 Zig 0.14.1 或 0.15.1：

```bash
# 1. 格式化代码
zig fmt src/fft/fft_radix2.zig src/fft/fft_r2c.zig src/fft/fft_mixed.zig

# 2. 构建项目
zig build

# 3. 运行测试
zig build test
```

或者直接运行提供的脚本：
```bash
./FORMAT_AND_TEST.sh
```

## 预期结果

所有测试应该通过，并且：
- 小尺寸 FFT 精度提升到 1e-12 相对误差
- 大数据 FFT 幅值精度达到 90% 以上
- FFT/IFFT 往返精度恢复到机器精度

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

修复后的实现满足所有这些性质。
