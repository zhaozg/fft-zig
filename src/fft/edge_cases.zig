//! Edge case and stress tests for FFT
//! 边界情况与压力测试

const std = @import("std");
const math = std.math;
const Complex = @import("types.zig").Complex;
const fft_module = @import("../fft.zig");
const fftInPlace = fft_module.fftInPlace;
const fft_radix2 = @import("fft_radix2.zig");
const fft_radix4 = @import("fft_radix4.zig");
const fft_mixed = @import("fft_mixed.zig");

const expectApproxEqRel = std.testing.expectApproxEqRel;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;
const expect = std.testing.expect;
const TEST_TOLERANCE = 1e-10;

// 测试零输入
test "Edge case: Zero input" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const sizes = [_]usize{ 4, 8, 16, 64 };
    for (sizes) |size| {
        var data = try allocator.alloc(Complex, size);
        defer allocator.free(data);

        for (0..size) |i| {
            data[i] = Complex{ .re = 0.0, .im = 0.0 };
        }

        try fftInPlace(allocator, data);

        // FFT of zeros should be zeros
        for (data) |v| {
            try expectApproxEqAbs(@as(f64, 0.0), v.re, TEST_TOLERANCE);
            try expectApproxEqAbs(@as(f64, 0.0), v.im, TEST_TOLERANCE);
        }
    }
}

// 测试极小值
test "Edge case: Very small values" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const size: usize = 64;
    var data = try allocator.alloc(Complex, size);
    defer allocator.free(data);

    const tiny_value: f64 = 1e-100;
    for (0..size) |i| {
        data[i] = Complex{ .re = tiny_value, .im = 0.0 };
    }

    try fftInPlace(allocator, data);

    // Should get DC component only
    const expected_dc = @as(f64, @floatFromInt(size)) * tiny_value;
    try expectApproxEqRel(expected_dc, data[0].re, 1e-6);
}

// 测试极大值
test "Edge case: Very large values" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const size: usize = 64;
    var data = try allocator.alloc(Complex, size);
    defer allocator.free(data);

    const large_value: f64 = 1e100;
    for (0..size) |i| {
        data[i] = Complex{ .re = large_value, .im = 0.0 };
    }

    try fftInPlace(allocator, data);

    // Should get DC component only
    const expected_dc = @as(f64, @floatFromInt(size)) * large_value;
    const rel_error = @abs(data[0].re - expected_dc) / expected_dc;
    try expect(rel_error < 1e-10);
}

// 测试复数输入
test "Edge case: Complex input" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const size: usize = 8;
    var data = try allocator.alloc(Complex, size);
    defer allocator.free(data);

    // 创建复数输入
    var prng = std.Random.DefaultPrng.init(777);
    const random = prng.random();
    for (0..size) |i| {
        data[i] = Complex{
            .re = random.float(f64) * 2.0 - 1.0,
            .im = random.float(f64) * 2.0 - 1.0,
        };
    }

    const original = try allocator.dupe(Complex, data);
    defer allocator.free(original);

    try fftInPlace(allocator, data);

    // 验证能量守恒 (Parseval)
    var time_energy: f64 = 0.0;
    for (original) |v| {
        time_energy += v.re * v.re + v.im * v.im;
    }

    var freq_energy: f64 = 0.0;
    for (data) |v| {
        freq_energy += v.re * v.re + v.im * v.im;
    }

    const expected_freq_energy = time_energy * @as(f64, @floatFromInt(size));
    try expectApproxEqRel(expected_freq_energy, freq_energy, 0.01);
}

// 测试交替符号输入
test "Edge case: Alternating signs" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const size: usize = 16;
    var data = try allocator.alloc(Complex, size);
    defer allocator.free(data);

    // [1, -1, 1, -1, ...]
    for (0..size) |i| {
        data[i] = Complex{
            .re = if (i % 2 == 0) 1.0 else -1.0,
            .im = 0.0,
        };
    }

    try fftInPlace(allocator, data);

    // 应该在 N/2 处有峰值
    const nyquist_bin = size / 2;
    var max_mag: f64 = 0.0;
    var max_bin: usize = 0;

    for (0..size) |k| {
        const mag = @sqrt(data[k].re * data[k].re + data[k].im * data[k].im);
        if (mag > max_mag) {
            max_mag = mag;
            max_bin = k;
        }
    }

    try expect(max_bin == nyquist_bin);
}

// 测试不同大小的非2的幂
test "Edge case: Non-power-of-2 sizes" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const test_sizes = [_]usize{ 3, 5, 6, 7, 9, 10, 12, 15 };

    for (test_sizes) |size| {
        var data = try allocator.alloc(Complex, size);
        defer allocator.free(data);

        // 生成简单测试信号
        for (0..size) |i| {
            data[i] = Complex{ .re = @as(f64, @floatFromInt(i + 1)), .im = 0.0 };
        }

        // 应该使用 mixed-radix 或 Bluestein
        try fftInPlace(allocator, data);

        // 基本正确性：DC分量检查
        var expected_dc: f64 = 0.0;
        for (0..size) |i| {
            expected_dc += @as(f64, @floatFromInt(i + 1));
        }
        try expectApproxEqRel(expected_dc, data[0].re, 0.01);
    }
}

// 测试精度：比较FFT与DFT
test "Edge case: Accuracy comparison FFT vs DFT" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const sizes = [_]usize{ 8, 16, 32 };

    for (sizes) |size| {
        var fft_data = try allocator.alloc(Complex, size);
        defer allocator.free(fft_data);
        var dft_data = try allocator.alloc(Complex, size);
        defer allocator.free(dft_data);

        // 生成测试信号
        var prng = std.Random.DefaultPrng.init(888);
        const random = prng.random();
        for (0..size) |i| {
            const val = Complex{
                .re = random.float(f64) * 2.0 - 1.0,
                .im = 0.0,
            };
            fft_data[i] = val;
            dft_data[i] = val;
        }

        // FFT
        try fftInPlace(allocator, fft_data);

        // DFT (直接实现)
        var temp = try allocator.alloc(Complex, size);
        defer allocator.free(temp);
        for (0..size) |k| {
            temp[k] = Complex{ .re = 0.0, .im = 0.0 };
            for (0..size) |n| {
                const angle = -2.0 * math.pi * @as(f64, @floatFromInt(k * n)) / @as(f64, @floatFromInt(size));
                const w = Complex{ .re = @cos(angle), .im = @sin(angle) };
                temp[k].re += dft_data[n].re * w.re - dft_data[n].im * w.im;
                temp[k].im += dft_data[n].re * w.im + dft_data[n].im * w.re;
            }
        }

        // 比较结果
        for (0..size) |k| {
            // Use absolute tolerance for very small values to avoid division by near-zero
            const abs_val_re = @abs(temp[k].re);
            const abs_val_im = @abs(temp[k].im);

            if (abs_val_re < 1e-12) {
                try expectApproxEqAbs(temp[k].re, fft_data[k].re, 1e-10);
            } else {
                try expectApproxEqRel(temp[k].re, fft_data[k].re, 1e-8);
            }

            if (abs_val_im < 1e-12) {
                try expectApproxEqAbs(temp[k].im, fft_data[k].im, 1e-10);
            } else {
                try expectApproxEqRel(temp[k].im, fft_data[k].im, 1e-8);
            }
        }
    }
}

// 测试多次调用的一致性
test "Edge case: Consistency across multiple calls" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const size: usize = 64;
    var original = try allocator.alloc(Complex, size);
    defer allocator.free(original);

    // 生成测试信号
    var prng = std.Random.DefaultPrng.init(333);
    const random = prng.random();
    for (0..size) |i| {
        original[i] = Complex{
            .re = random.float(f64) * 2.0 - 1.0,
            .im = 0.0,
        };
    }

    var result1 = try allocator.dupe(Complex, original);
    defer allocator.free(result1);
    var result2 = try allocator.dupe(Complex, original);
    defer allocator.free(result2);
    var result3 = try allocator.dupe(Complex, original);
    defer allocator.free(result3);

    // 多次调用应该得到相同结果
    try fftInPlace(allocator, result1);
    try fftInPlace(allocator, result2);
    try fftInPlace(allocator, result3);

    // Explicitly use mutable references to satisfy compiler
    _ = &result1;
    _ = &result2;
    _ = &result3;

    for (0..size) |i| {
        try expectApproxEqRel(result1[i].re, result2[i].re, TEST_TOLERANCE);
        try expectApproxEqRel(result1[i].im, result2[i].im, TEST_TOLERANCE);
        try expectApproxEqRel(result2[i].re, result3[i].re, TEST_TOLERANCE);
        try expectApproxEqRel(result2[i].im, result3[i].im, TEST_TOLERANCE);
    }
}

// 测试边界条件：size = 1
test "Edge case: Size 1" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var data = [_]Complex{Complex{ .re = 3.14159, .im = 2.71828 }};
    const original = data[0];

    try fftInPlace(allocator, &data);

    // Size 1 FFT should be identity
    try expectApproxEqRel(original.re, data[0].re, TEST_TOLERANCE);
    try expectApproxEqRel(original.im, data[0].im, TEST_TOLERANCE);
}

// 测试特定模式：脉冲序列
test "Edge case: Impulse train" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const size: usize = 16;
    var data = try allocator.alloc(Complex, size);
    defer allocator.free(data);

    // 每4个样本一个脉冲
    for (0..size) |i| {
        data[i] = Complex{
            .re = if (i % 4 == 0) 1.0 else 0.0,
            .im = 0.0,
        };
    }

    try fftInPlace(allocator, data);

    // 验证周期性 - 脉冲序列在频域也应该是周期性的
    // 对于周期为4的脉冲序列，频域应该在k=0, 4, 8, 12处有值
    const period = 4;
    for (0..size) |k| {
        const mag = @sqrt(data[k].re * data[k].re + data[k].im * data[k].im);
        if (k % period == 0) {
            // 应该有非零值
            try expect(mag > 0.1);
        }
    }
}
