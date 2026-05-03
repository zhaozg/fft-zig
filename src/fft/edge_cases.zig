//! Edge case and stress tests for FFT
//! 边界情况与压力测试

const std = @import("std");
const math = std.math;
const fft_module = @import("../fft.zig");
const fftInPlace = fft_module.fftInPlace;

const expectApproxEqRel = std.testing.expectApproxEqRel;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;
const expect = std.testing.expect;

fn tolerance(comptime T: type) T {
    return if (T == f32) @as(T, 1e-4) else @as(T, 1e-10);
}

// 测试零输入
fn testZeroInputGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tol = tolerance(T);

    const sizes = [_]usize{ 4, 8, 16, 64 };
    for (sizes) |size| {
        var data = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(data);

        for (0..size) |i| {
            data[i] = std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) };
        }

        try fftInPlace(T, allocator, data);

        for (data) |v| {
            try expectApproxEqAbs(@as(T, 0.0), v.re, tol);
            try expectApproxEqAbs(@as(T, 0.0), v.im, tol);
        }
    }
}

test "Edge case: Zero input f32" {
    try testZeroInputGeneric(f32);
}
test "Edge case: Zero input f64" {
    try testZeroInputGeneric(f64);
}

// 测试极小值
fn testTinyValuesGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;

    const size: usize = 64;
    var data = try allocator.alloc(std.math.Complex(T), size);
    defer allocator.free(data);

    const tiny_value: T = @as(T, 1e-30);
    for (0..size) |i| {
        data[i] = std.math.Complex(T){ .re = tiny_value, .im = @as(T, 0.0) };
    }

    try fftInPlace(T, allocator, data);

    const expected_dc = @as(T, @floatFromInt(size)) * tiny_value;
    try expectApproxEqRel(expected_dc, data[0].re, @as(T, 1e-5));
}

test "Edge case: Very small values f32" {
    try testTinyValuesGeneric(f32);
}
test "Edge case: Very small values f64" {
    try testTinyValuesGeneric(f64);
}

// 测试复数输入
fn testComplexInputGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;

    const size: usize = 8;
    var data = try allocator.alloc(std.math.Complex(T), size);
    defer allocator.free(data);

    var prng = std.Random.DefaultPrng.init(777);
    const random = prng.random();
    for (0..size) |i| {
        data[i] = std.math.Complex(T){
            .re = random.float(T) * @as(T, 2.0) - @as(T, 1.0),
            .im = random.float(T) * @as(T, 2.0) - @as(T, 1.0),
        };
    }

    const original = try allocator.dupe(std.math.Complex(T), data);
    defer allocator.free(original);

    try fftInPlace(T, allocator, data);

    var time_energy: T = @as(T, 0.0);
    for (original) |v| {
        time_energy += v.re * v.re + v.im * v.im;
    }

    var freq_energy: T = @as(T, 0.0);
    for (data) |v| {
        freq_energy += v.re * v.re + v.im * v.im;
    }

    const expected_freq_energy = time_energy * @as(T, @floatFromInt(size));
    try expectApproxEqRel(expected_freq_energy, freq_energy, @as(T, 0.01));
}

test "Edge case: Complex input f32" {
    try testComplexInputGeneric(f32);
}
test "Edge case: Complex input f64" {
    try testComplexInputGeneric(f64);
}

// 测试交替符号输入
fn testAlternatingSignsGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;

    const size: usize = 16;
    var data = try allocator.alloc(std.math.Complex(T), size);
    defer allocator.free(data);

    for (0..size) |i| {
        data[i] = std.math.Complex(T){
            .re = if (i % 2 == 0) @as(T, 1.0) else @as(T, -1.0),
            .im = @as(T, 0.0),
        };
    }

    try fftInPlace(T, allocator, data);

    const nyquist_bin = size / 2;
    var max_mag: T = @as(T, 0.0);
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

test "Edge case: Alternating signs f32" {
    try testAlternatingSignsGeneric(f32);
}
test "Edge case: Alternating signs f64" {
    try testAlternatingSignsGeneric(f64);
}

// 测试不同大小的非2的幂
fn testNonPowerOfTwoGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;

    const test_sizes = [_]usize{ 3, 5, 6, 7, 9, 10, 12, 15 };

    for (test_sizes) |size| {
        var data = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(data);

        for (0..size) |i| {
            data[i] = std.math.Complex(T){ .re = @as(T, @floatFromInt(i + 1)), .im = @as(T, 0.0) };
        }

        try fftInPlace(T, allocator, data);

        var expected_dc: T = @as(T, 0.0);
        for (0..size) |i| {
            expected_dc += @as(T, @floatFromInt(i + 1));
        }
        try expectApproxEqRel(expected_dc, data[0].re, @as(T, 0.01));
    }
}

test "Edge case: Non-power-of-2 sizes f32" {
    try testNonPowerOfTwoGeneric(f32);
}
test "Edge case: Non-power-of-2 sizes f64" {
    try testNonPowerOfTwoGeneric(f64);
}

// 测试精度：比较FFT与DFT
fn testAccuracyGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tol = tolerance(T);

    const sizes = [_]usize{ 8, 16, 32 };

    for (sizes) |size| {
        var fft_data = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(fft_data);
        var dft_data = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(dft_data);

        var prng = std.Random.DefaultPrng.init(888);
        const random = prng.random();
        for (0..size) |i| {
            const val = std.math.Complex(T){
                .re = random.float(T) * @as(T, 2.0) - @as(T, 1.0),
                .im = @as(T, 0.0),
            };
            fft_data[i] = val;
            dft_data[i] = val;
        }

        try fftInPlace(T, allocator, fft_data);

        var temp = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(temp);
        for (0..size) |k| {
            temp[k] = std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) };
            for (0..size) |n| {
                const angle = -2.0 * math.pi * @as(T, @floatFromInt(k * n)) / @as(T, @floatFromInt(size));
                const w = std.math.Complex(T){ .re = @cos(angle), .im = @sin(angle) };
                temp[k].re += dft_data[n].re * w.re - dft_data[n].im * w.im;
                temp[k].im += dft_data[n].re * w.im + dft_data[n].im * w.re;
            }
        }

        for (0..size) |k| {
            try expectApproxEqAbs(temp[k].re, fft_data[k].re, tol);
            try expectApproxEqAbs(temp[k].im, fft_data[k].im, tol);
        }
    }
}

test "Edge case: Accuracy comparison FFT vs DFT f32" {
    try testAccuracyGeneric(f32);
}
test "Edge case: Accuracy comparison FFT vs DFT f64" {
    try testAccuracyGeneric(f64);
}

// 测试多次调用的一致性
fn testConsistencyGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tol = tolerance(T);

    const size: usize = 64;
    var original = try allocator.alloc(std.math.Complex(T), size);
    defer allocator.free(original);

    var prng = std.Random.DefaultPrng.init(333);
    const random = prng.random();
    for (0..size) |i| {
        original[i] = std.math.Complex(T){
            .re = random.float(T) * @as(T, 2.0) - @as(T, 1.0),
            .im = @as(T, 0.0),
        };
    }

    const result1 = try allocator.dupe(std.math.Complex(T), original);
    defer allocator.free(result1);
    const result2 = try allocator.dupe(std.math.Complex(T), original);
    defer allocator.free(result2);

    try fftInPlace(T, allocator, result1);
    try fftInPlace(T, allocator, result2);

    for (0..size) |i| {
        try expectApproxEqRel(result1[i].re, result2[i].re, tol);
        try expectApproxEqRel(result1[i].im, result2[i].im, tol);
    }
}

test "Edge case: Consistency across multiple calls f32" {
    try testConsistencyGeneric(f32);
}
test "Edge case: Consistency across multiple calls f64" {
    try testConsistencyGeneric(f64);
}

// 测试边界条件：size = 1
fn testSizeOneGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tol = tolerance(T);

    var data = [_]std.math.Complex(T){std.math.Complex(T){ .re = @as(T, 3.14159), .im = @as(T, 2.71828) }};
    const original = data[0];

    try fftInPlace(T, allocator, &data);

    try expectApproxEqRel(original.re, data[0].re, tol);
    try expectApproxEqRel(original.im, data[0].im, tol);
}

test "Edge case: Size 1 f32" {
    try testSizeOneGeneric(f32);
}
test "Edge case: Size 1 f64" {
    try testSizeOneGeneric(f64);
}

// 测试特定模式：脉冲序列
fn testImpulseTrainGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;

    const size: usize = 16;
    var data = try allocator.alloc(std.math.Complex(T), size);
    defer allocator.free(data);

    for (0..size) |i| {
        data[i] = std.math.Complex(T){
            .re = if (i % 4 == 0) @as(T, 1.0) else @as(T, 0.0),
            .im = @as(T, 0.0),
        };
    }

    try fftInPlace(T, allocator, data);

    const period = 4;
    for (0..size) |k| {
        const mag = @sqrt(data[k].re * data[k].re + data[k].im * data[k].im);
        if (k % period == 0) {
            try expect(mag > @as(T, 0.1));
        }
    }
}

test "Edge case: Impulse train f32" {
    try testImpulseTrainGeneric(f32);
}
test "Edge case: Impulse train f64" {
    try testImpulseTrainGeneric(f64);
}
