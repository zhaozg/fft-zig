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

## 技术说明

FFT 本身只是一个数学算法，它不关心信号的物理意义。
当我们对一个时域信号 `x(n)` 进行 FFT 得到频域结果 `X(k)` 时，`X(k)` 的幅度和功率
与原始信号的采样点数 `N` 直接相关。为了使得频域结果（如幅度、功率）具有明确的物
理意义，并且与采样点数 `N` 无关，我们必须进行归一化。

### 1. FFT 的数学定义与归一化因子的关系

理解归一化的关键在于回顾 DFT（离散傅里叶变换）的数学定义。FFT 只是 DFT 的一种快
速算法，其数学基础是相同的。

**DFT 正变换（分析方程）：**
\[
X(k) = \sum_{n=0}^{N-1} x(n) \cdot e^{-j \frac{2\pi}{N}kn}
\]
**DFT 逆变换（合成方程）：**
\[
x(n) = \frac{1}{N} \sum_{k=0}^{N-1} X(k) \cdot e^{j \frac{2\pi}{N}kn}
\]

请注意，**归一化因子（通常是 `1/N`）可以放在正变换或逆变换中，也可以拆分开
（如 `1/sqrt(N)` 放在两边）**。不同的软件库和教科书有不同的约定，这是混淆的主要
来源。

### 2. 常见的归一化方法及其应用场景

下面我们以最常见的 **`1/N` 归一化** 为例，详细说明几种方法。

假设我们有一个实数信号，例如一个幅度为 `A`，频率为 `f` 的单频余弦波：
\[
x(n) = A \cdot \cos(2\pi f n / N)
\]
它的 FFT 结果在正负频率处各有一个峰。

#### 方法一：幅度归一化 - 用于分析信号幅度

这是**最常用**的方法，目的是让频域谱线的幅度直接对应原始信号中该频率分量的**实际幅度**。

*   **归一化因子：** `1 / N`
*   **应用位置：** 对 FFT 计算出的**全部**复数结果 `X(k)` 进行缩放。
    \[
    X\_{\text{normalized}}(k) = \frac{X(k)}{N}
    \]
*   **结果解释：**
    *   对于**复数信号**的单边频谱，峰值处的 `|X(k)|` 就是该频率分量的幅度 `A`。
    *   对于**实数信号**，其能量会分布在正负频率两个点上。因此：
        *   **双边频谱：** 在 `k=f` 和 `k=N-f` 处的两个峰，每个峰的幅度约为 `A/2`。归一化后，每个峰的 `|X(k)|` 约为 `A/2`。
        *   **单边频谱：** 将正频率部分的幅度乘以 `2`（负频率部分被丢弃），这样在 `k=f` 处的峰值幅度就恢复为 `A`。

*   **应用场景：**
    *   分析振动信号的位移、速度、加速度幅度。
    *   分析交流电压/电流的幅度。
    *   任何需要知道信号中各频率分量具体幅度的场合。

**示例：**
一个幅度为 0.8，频率为 50 Hz 的余弦波，做 `N=1024` 的 FFT。
- 未经归一化的 FFT 在 50Hz 和 (1024-50)Hz 处的峰值幅度约为 `~410`。
- 经过 `1/N` 归一化后，这两个峰的幅度变为 `~410 / 1024 ≈ 0.4`。
- 绘制单边频谱时，再将 50Hz 处的幅度乘以 2，得到 `0.8`，正好等于原始信号的幅度。

#### 方法二：功率归一化 - 用于分析信号功率

此方法用于还原信号的**功率**或**能量**。

*   **归一化因子：** `1 / N` 用于线性功率，`1 / N^2` 用于功率谱密度。
*   **常见形式：**
    1.  **功率谱：** 先计算 FFT 结果的平方（得到功率），然后除以 `N`。
        \[
        P(k) = \frac{|X(k)|^2}{N}
        \]
        或者，为了得到单边功率谱，在正频率部分还需要再乘以 `2`（负频率部分的功率被加到此）。
        \[
        P\_{\text{single-sided}}(k) = \frac{2 \cdot |X(k)|^2}{N} \quad (k=1,2,...,N/2-1)
        \]
        这样，`P(k)` 的数值就代表了该频率分量的功率。

    2.  **功率谱密度：** 为了比较不同采样率下的功率，需要除以频率分辨率 `Δf = Fs / N`，其中 `Fs` 是采样频率。这使 PSD 的单位是 `V²/Hz`。
        \[
        PSD(k) = \frac{|X(k)|^2}{N \cdot \Delta f} = \frac{|X(k)|^2 \cdot N}{F_s}
        \]
        （同样，单边 PSD 需要在正频率部分乘以 `2`）。

*   **应用场景：**
    *   计算信号的总功率。
    *   在通信系统中分析信噪比。
    *   随机振动分析。

#### 方法三：能量归一化 - 用于分析暂态信号能量

根据帕塞瓦尔定理，时域的总能量等于频域的总能量。

*   **定理：** \(\sum_{n=0}^{N-1} |x(n)|^2 = \frac{1}{N} \sum_{k=0}^{N-1} |X(k)|^2\)
*   **归一化方法：** 为了满足这个定理，通常的作法是在 FFT 之后除以 `N`，然后再计算平方。
    \[
    E = \sum |x(n)|^2 = \sum \frac{|X(k)|^2}{N}
    \]
    可以看到，这与功率归一化中的功率谱公式是一致的。所以能量归一化可以看作是功率归一化在特定上下文下的另一种解释。

*   **应用场景：**
    *   分析脉冲、冲击等暂态信号的总能量。

#### 方法四：对称归一化 (`1/sqrt(N)`)

这是一种在数学上很优雅的方法，它使正变换和逆变换在形式上完全对称。

*   **正变换：** \(X(k) = \frac{1}{\sqrt{N}} \sum_{n=0}^{N-1} x(n) \cdot e^{-j \frac{2\pi}{N}kn}\)
*   **逆变换：** \(x(n) = \frac{1}{\sqrt{N}} \sum_{k=0}^{N-1} X(k) \cdot e^{j \frac{2\pi}{N}kn}\)

*   **特点：**
    *   它保持了时域和频域的“范数”（总能量）不变。
    *   但在工程实际中不如 `1/N` 直观，因为峰值幅度不再直接等于 `A`。
*   **应用场景：**
    *   更多见于纯数学和物理领域的推导。

### 不同软件库的默认行为

了解你使用的工具库的默认行为至关重要：

*   **MATLAB / Octave:**
    *   `fft(x)` 函数**不做任何归一化**。你需要手动除以 `N` 来进行幅度归一化。
    *   `ifft(X)` 函数会自动除以 `N`。所以 `ifft(fft(x))` 会得到原始信号 `x`。

*   **Python (NumPy/SciPy):**
    *   `numpy.fft.fft(x)` 和 MATLAB 一样，**不做任何归一化**。
    *   `numpy.fft.ifft(X)` 会自动除以 `N`。

*   **Python (SciPy.signal):**
    *   `scipy.signal.spectrogram` 等函数通常提供 `scaling` 参数（如 `'density'` 或 `'spectrum'`）让你选择功率谱密度或功率谱归一化。

*  **fft-zig**
    *  **不做任何归一化**。

### 验证方法

可以用以下方法验证 FFT 正确性：

1. **Parseval 定理**: Σ|x[n]|² = (1/N)·Σ|X[k]|²
2. **单位冲激**: FFT([1,0,0,...,0]) = [1,1,1,...,1]
3. **单频正弦**: FFT(sin(2πk₀n/N)) 在 k=k₀ 处有峰值 N/2

## 贡献 / Contributing

- 代码使用 `zig fmt` 格式化 / Code formatted with `zig fmt`

## 许可 / License

See repository license
