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
const fft_r2c = @import("fft/fft_r2c.zig");
const fft_utils = @import("fft/utils.zig");
const fft_inplace = @import("fft/fft_inplace.zig");

// 向后兼容的类型别名
pub const Complex = fft_types.ComplexF64;
pub const TwiddleFactorTable = fft_twiddle.TwiddleFactorTable;
const isPowerOfTwo = fft_utils.isPowerOfTwo;
const isPowerOfFour = fft_utils.isPowerOfFour;
pub const fftInPlace = fft_inplace.fftInPlace;

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
