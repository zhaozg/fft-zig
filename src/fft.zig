//! High-performance FFT (Fast Fourier Transform) implementation
//! Optimized for performance with SIMD vectorization and parallel processing
//!
//! This module provides a comprehensive FFT implementation targeting GSL-level performance
//! Features:
//! - SIMD-optimized radix-2, radix-4, and mixed-radix algorithms
//! - Compile-time optimized twiddle factor tables
//! - Parallel processing for large datasets
//! - Automatic algorithm selection based on input size
//! - Support for both power-of-2 and arbitrary-length transforms

const std = @import("std");
const math = std.math;
const fft_types = @import("fft/types.zig");
const fft_twiddle = @import("fft/twiddle.zig");
const fft_radix2 = @import("fft/fft_radix2.zig");
const fft_radix4 = @import("fft/fft_radix4.zig");
const fft_mixed = @import("fft/fft_mixed.zig");
const fft_parallel = @import("fft/fft_parallel.zig");
const fft_r2c = @import("fft/fft_r2c.zig");
const fft_utils = @import("fft/utils.zig");

// 向后兼容的类型别名
pub const Complex = fft_types.ComplexF64;
pub const TwiddleFactorTable = fft_twiddle.TwiddleFactorTable;
const isPowerOfTwo = fft_utils.isPowerOfTwo;
const isPowerOfFour = fft_utils.isPowerOfFour;

// Performance thresholds for algorithm selection
pub const PARALLEL_THRESHOLD = 16384;
pub const HUGE_DATA_THRESHOLD = 1000000; // 1M threshold for huge data optimizations
pub const SIMD_THRESHOLD = 64;
pub const RADIX4_THRESHOLD = 256;
pub const SMALL_FFT_THRESHOLD = 256;

/// 泛型 FFT 接口函数（实数到复数）
/// T: 浮点类型（f32 或 f64）
/// allocator: 内存分配器
/// input: 实数输入数据
/// output: 复数输出数据
pub fn fft(comptime T: type, allocator: std.mem.Allocator, input: []const T, output: []std.math.Complex(T)) !void {
    // 分配临时缓冲区
    const tmp_out = try allocator.alloc(T, 2 * output.len);
    defer allocator.free(tmp_out);

    // 幅值谱可选
    const tmp_mag = try allocator.alloc(T, output.len);
    defer allocator.free(tmp_mag);

    // 调用高性能 FFT
    try fft_r2c.fftR2C(T, allocator, input, tmp_out, tmp_mag);

    // 转换为 Complex 输出
    output[0] = std.math.Complex(T){ .re = tmp_out[0], .im = @as(T, 0.0) };
    for (1..output.len) |k| {
        output[k] = std.math.Complex(T){ .re = tmp_out[2 * k], .im = tmp_out[2 * k + 1] };
    }
}

/// 泛型就地 FFT 函数
pub fn fftInPlace(comptime T: type, allocator: std.mem.Allocator, data: []std.math.Complex(T)) error{ InvalidSize, OutOfMemory }!void {
    // 自动选择算法，委托给各模块
    const n = data.len;
    if (n >= HUGE_DATA_THRESHOLD) {
        try fft_parallel.fftHugeDataParallel(T, allocator, data);
    } else if (n >= PARALLEL_THRESHOLD and fft_utils.isPowerOfTwo(n)) {
        try fft_parallel.fftParallelSIMD(T, allocator, data);
    } else if (n >= RADIX4_THRESHOLD and fft_utils.isPowerOfFour(n)) {
        try fft_radix4.fftRadix4SIMD(T, data);
    } else if (n >= SIMD_THRESHOLD and fft_utils.isPowerOfTwo(n)) {
        try fft_radix2.fftRadix2SIMD(T, data);
    } else if (fft_utils.isPowerOfFour(n)) {
        try fft_radix4.fftRadix4(T, data);
    } else if (fft_utils.isPowerOfTwo(n)) {
        try fft_radix2.fftRadix2(T, data);
    } else if (n % 4 == 0 and fft_utils.isPowerOfTwo(n / 4)) {
        try fft_radix4.fftRadix4(T, data);
    } else {
        try fft_mixed.fftMixedRadix(T, allocator, data);
    }
}

/// Direct DFT implementation for testing and verification
pub fn dft(comptime T: type, input: []const std.math.Complex(T), output: []std.math.Complex(T)) void {
    const n = input.len;

    for (0..n) |k| {
        output[k] = std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) };

        for (0..n) |j| {
            const angle = -2.0 * math.pi * @as(T, @floatFromInt(k)) * @as(T, @floatFromInt(j)) / @as(T, @floatFromInt(n));
            const w = std.math.Complex(T){
                .re = math.cos(angle),
                .im = math.sin(angle),
            };

            output[k] = output[k].add(input[j].mul(w));
        }
    }
}

// 向后兼容的便捷别名
pub const fft_f64 = fft(f64);
pub const fft_f32 = fft(f32);
pub const fftInPlace_f64 = fftInPlace(f64);
pub const fftInPlace_f32 = fftInPlace(f32);

// 导出 fftR2C 供外部使用
pub const fftR2C = fft_r2c.fftR2C;

const expectApproxEqRel = std.testing.expectApproxEqRel;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;
const expect = std.testing.expect;

fn testFFTPerformanceGeneric(comptime T: type) !void {
    const io = std.testing.io;
    const allocator = std.testing.allocator;

    const test_sizes = [_]usize{ 256, 1024, 4096 };
    for (test_sizes) |size| {
        std.debug.print("Benchmarking size {d}...\n", .{size});

        var input = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(input);

        for (0..size) |i| {
            const t = @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(size));
            input[i] = std.math.Complex(T){ .re = @sin(2.0 * std.math.pi * @as(T, 7.0) * t) + @as(T, 0.3) * @cos(2.0 * std.math.pi * @as(T, 23.0) * t), .im = @as(T, 0.0) };
        }

        const data = try allocator.dupe(std.math.Complex(T), input);
        defer allocator.free(data);

        const clock = std.Io.Clock.awake;
        const start_time = std.Io.Clock.now(clock, io).toNanoseconds();
        try fftInPlace(T, allocator, data);
        const end_time = std.Io.Clock.now(clock, io).toNanoseconds();

        const elapsed_s = @as(f128, @floatFromInt(@as(i96, @intCast(end_time - start_time)))) / std.time.us_per_s;
        const throughput = (@as(f128, @floatFromInt(size)) / (elapsed_s) / 1e6);

        std.debug.print("  Size {d}: {d:.2}ms, {d:.1} MSamples/s\n", .{ size, elapsed_s, throughput });

        // Validate correctness - find dominant frequency
        var max_magnitude: T = @as(T, 0.0);
        var peak_freq: usize = 0;
        for (1..size / 2) |i| {
            const magnitude = @sqrt(data[i].re * data[i].re + data[i].im * data[i].im);
            if (magnitude > max_magnitude) {
                max_magnitude = magnitude;
                peak_freq = i;
            }
        }
        const expected_freq = (7 * size) / size;
        try expect(peak_freq >= expected_freq - 2 and peak_freq <= expected_freq + 2);
    }
}

test "FFT performance and large data f32" { try testFFTPerformanceGeneric(f32); }
test "FFT performance and large data f64" { try testFFTPerformanceGeneric(f64); }

test {
    _ = @import("fft/base.zig");
    _ = @import("fft/fft_r2c.zig");
    _ = @import("fft/ifft.zig");
    _ = @import("fft/utils.zig");
    _ = @import("fft/fft_mixed.zig");
    _ = @import("fft/fft_radix2.zig");
    _ = @import("fft/twiddle.zig");
    _ = @import("fft/fft_parallel.zig");
    _ = @import("fft/fft_radix4.zig");
    _ = @import("fft/types.zig");
    _ = @import("fft/validation.zig");
    _ = @import("fft/performance.zig");
    _ = @import("fft/edge_cases.zig");
}
