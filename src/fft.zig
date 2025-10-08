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
const builtin = @import("builtin");
const fft_types = @import("fft/types.zig");
const fft_twiddle = @import("fft/twiddle.zig");
const fft_radix2 = @import("fft/fft_radix2.zig");
const fft_radix4 = @import("fft/fft_radix4.zig");
const fft_mixed = @import("fft/fft_mixed.zig");
const fft_parallel = @import("fft/fft_parallel.zig");
const fft_ifft = @import("fft/ifft.zig");
const fft_r2c = @import("fft/fft_r2c.zig");
const fft_utils = @import("fft/utils.zig");

pub const fftR2C = fft_r2c.fftR2C;

// 类型定义和 SIMD 支持已迁移到 src/fft/types.zig
const VectorF64 = fft_types.VectorF64;
pub const Complex = fft_types.Complex;
pub const TwiddleFactorTable = fft_twiddle.TwiddleFactorTable;
const isPowerOfTwo = fft_utils.isPowerOfTwo;
const isPowerOfFour = fft_utils.isPowerOfFour;
const nextPowerOfTwo = fft_utils.nextPowerOfTwo;

// Performance thresholds for algorithm selection
pub const PARALLEL_THRESHOLD = 16384;
pub const HUGE_DATA_THRESHOLD = 1000000; // 1M threshold for huge data optimizations
pub const SIMD_THRESHOLD = 64;
pub const RADIX4_THRESHOLD = 256;
pub const SMALL_FFT_THRESHOLD = 256;

// Twiddle 表和相关工具已迁移到 src/fft/twiddle.zig

/// FFT 接口函数，具体实现已迁移到 src/fft 目录
pub fn fft(allocator: std.mem.Allocator, input: []const f64, output: []Complex) !void {
    // 分配临时缓冲区
    const tmp_out = try allocator.alloc(f64, 2 * output.len);
    defer allocator.free(tmp_out);

    // 幅值谱可选
    const tmp_mag = try allocator.alloc(f64, output.len);
    defer allocator.free(tmp_mag);

    // 调用高性能 FFT
    try fft_r2c.fftR2C(allocator, input, tmp_out, tmp_mag);

    // 转换为 Complex 输出
    output[0] = Complex{ .re = tmp_out[0], .im = 0.0 };
    for (1..output.len) |k| {
        output[k] = Complex{ .re = tmp_out[2 * k], .im = tmp_out[2 * k + 1] };
    }
}

pub fn fftInPlace(allocator: std.mem.Allocator, data: []Complex) error{ InvalidSize, OutOfMemory }!void {
    // 自动选择算法，委托给各模块
    const n = data.len;
    if (n >= HUGE_DATA_THRESHOLD) {
        try fft_parallel.fftHugeDataParallel(allocator, data);
    } else if (n >= PARALLEL_THRESHOLD and fft_utils.isPowerOfTwo(n)) {
        try fft_parallel.fftParallelSIMD(allocator, data);
    } else if (n >= RADIX4_THRESHOLD and fft_utils.isPowerOfFour(n)) {
        try fft_radix4.fftRadix4SIMD(data);
    } else if (n >= SIMD_THRESHOLD and fft_utils.isPowerOfTwo(n)) {
        try fft_radix2.fftRadix2SIMD(data);
    } else if (fft_utils.isPowerOfFour(n)) {
        try fft_radix4.fftRadix4(data);
    } else if (fft_utils.isPowerOfTwo(n)) {
        try fft_radix2.fftRadix2(data);
    } else if (n % 4 == 0 and fft_utils.isPowerOfTwo(n / 4)) {
        try fft_radix4.fftRadix4(data);
    } else {
        try fft_mixed.fftMixedRadix(allocator, data);
    }
}

/// Direct DFT implementation for testing and verification
pub fn dft(input: []const Complex, output: []Complex) void {
    const n = input.len;

    for (0..n) |k| {
        output[k] = Complex{ .re = 0.0, .im = 0.0 };

        for (0..n) |j| {
            const angle = -2.0 * math.pi * @as(f64, @floatFromInt(k)) * @as(f64, @floatFromInt(j)) / @as(f64, @floatFromInt(n));
            const w = Complex{
                .re = math.cos(angle),
                .im = math.sin(angle),
            };

            output[k] = output[k].add(input[j].mul(w));
        }
    }
}

const expectApproxEqRel = std.testing.expectApproxEqRel;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;
const expect = std.testing.expect;
const TEST_TOLERANCE = 1e-12;

test "FFT performance and large data" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const test_sizes = [_]usize{ 256, 1024, 4096 };
    for (test_sizes) |size| {
        std.debug.print("Benchmarking size {d}...\n", .{size});

        var input = try allocator.alloc(Complex, size);
        defer allocator.free(input);

        for (0..size) |i| {
            const t = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(size));
            input[i] = Complex{ .re = @sin(2.0 * std.math.pi * 7.0 * t) + 0.3 * @cos(2.0 * std.math.pi * 23.0 * t), .im = 0.0 };
        }

        const data = try allocator.dupe(Complex, input);
        defer allocator.free(data);

        const start_time = std.time.nanoTimestamp();
        try fftInPlace(allocator, data);
        const end_time = std.time.nanoTimestamp();

        const elapsed_ms = @as(f64, @floatFromInt(@as(u64, @intCast(end_time - start_time)))) / 1e6;
        const throughput = (@as(f64, @floatFromInt(size)) / (elapsed_ms / 1000.0)) / 1e6;

        std.debug.print("  Size {d}: {d:.2}ms, {d:.1} MSamples/s\n", .{ size, elapsed_ms, throughput });

        // Validate correctness - find dominant frequency
        var max_magnitude: f64 = 0.0;
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
