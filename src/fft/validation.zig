//! Comprehensive FFT validation tests
//! 全面的FFT正确性验证测试

const std = @import("std");
const math = std.math;
const fftInPlace = @import("base.zig").fftInPlaceBase;
const ifft_module = @import("ifft.zig");
const ifftInPlace = ifft_module.ifftInPlace;

const expectApproxEqRel = std.testing.expectApproxEqRel;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;
const expect = std.testing.expect;

fn tolerance(comptime T: type) T {
    return if (T == f32) @as(T, 1e-3) else @as(T, 1e-10);
}

// 测试1: 验证单位冲激响应
fn testUnitImpulseGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tol = tolerance(T);

    const sizes = [_]usize{ 4, 8, 16, 64, 256, 1024 };
    for (sizes) |size| {
        var data = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(data);

        data[0] = std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) };
        for (1..size) |i| {
            data[i] = std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) };
        }

        try fftInPlace(T, allocator, data);

        for (data) |v| {
            try expectApproxEqRel(@as(T, 1.0), v.re, tol);
            try expectApproxEqAbs(@as(T, 0.0), v.im, tol);
        }
    }
}

test "FFT validation: Unit impulse response f32" { try testUnitImpulseGeneric(f32); }
test "FFT validation: Unit impulse response f64" { try testUnitImpulseGeneric(f64); }

// 测试2: 验证直流分量
fn testDCComponentGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tol = tolerance(T);

    const sizes = [_]usize{ 4, 8, 16, 64, 256 };
    for (sizes) |size| {
        var data = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(data);

        const dc_value: T = @as(T, 3.14159);
        for (0..size) |i| {
            data[i] = std.math.Complex(T){ .re = dc_value, .im = @as(T, 0.0) };
        }

        try fftInPlace(T, allocator, data);

        const expected_dc = @as(T, @floatFromInt(size)) * dc_value;
        try expectApproxEqRel(expected_dc, data[0].re, tol);
        try expectApproxEqAbs(@as(T, 0.0), data[0].im, tol);

        for (1..size) |k| {
            const mag = @sqrt(data[k].re * data[k].re + data[k].im * data[k].im);
            try expect(mag < @as(T, 1e-6));
        }
    }
}

test "FFT validation: DC component f32" { try testDCComponentGeneric(f32); }
test "FFT validation: DC component f64" { try testDCComponentGeneric(f64); }

// 测试3: 验证单频正弦波
fn testSingleFreqGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;

    const size: usize = 64;
    const freq: usize = 5;

    var data = try allocator.alloc(std.math.Complex(T), size);
    defer allocator.free(data);

    for (0..size) |n| {
        const angle = 2.0 * math.pi * @as(T, @floatFromInt(freq * n)) / @as(T, @floatFromInt(size));
        data[n] = std.math.Complex(T){ .re = @sin(angle), .im = @as(T, 0.0) };
    }

    try fftInPlace(T, allocator, data);

    var max_mag: T = @as(T, 0.0);
    var max_bin: usize = 0;
    for (0..size / 2) |k| {
        const mag = @sqrt(data[k].re * data[k].re + data[k].im * data[k].im);
        if (mag > max_mag) {
            max_mag = mag;
            max_bin = k;
        }
    }

    try expect(max_bin == freq);

    const expected_peak = @as(T, @floatFromInt(size)) / @as(T, 2.0);
    try expectApproxEqRel(expected_peak, max_mag, @as(T, 0.01));
}

test "FFT validation: Single frequency sine wave f32" { try testSingleFreqGeneric(f32); }
test "FFT validation: Single frequency sine wave f64" { try testSingleFreqGeneric(f64); }

// 测试4: 验证Parseval定理
fn testParsevalGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;

    const sizes = [_]usize{ 4, 8, 16, 64, 256, 1024 };
    var prng = std.Random.DefaultPrng.init(12345);
    const random = prng.random();

    for (sizes) |size| {
        var data = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(data);

        var time_energy: T = @as(T, 0.0);
        for (0..size) |i| {
            const val = random.float(T) * @as(T, 2.0) - @as(T, 1.0);
            data[i] = std.math.Complex(T){ .re = val, .im = @as(T, 0.0) };
            time_energy += val * val;
        }

        try fftInPlace(T, allocator, data);

        var freq_energy: T = @as(T, 0.0);
        for (data) |v| {
            freq_energy += v.re * v.re + v.im * v.im;
        }

        const expected_freq_energy = time_energy * @as(T, @floatFromInt(size));
        try expectApproxEqRel(expected_freq_energy, freq_energy, @as(T, 0.01));
    }
}

test "FFT validation: Parseval's theorem f32" { try testParsevalGeneric(f32); }
test "FFT validation: Parseval's theorem f64" { try testParsevalGeneric(f64); }

// 测试5: 验证FFT-IFFT往返
fn testRoundTripGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tol = tolerance(T);

    const sizes = [_]usize{ 4, 8, 16, 64, 256, 1024 };
    var prng = std.Random.DefaultPrng.init(54321);
    const random = prng.random();

    for (sizes) |size| {
        var original = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(original);
        var data = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(data);

        for (0..size) |i| {
            original[i] = std.math.Complex(T){
                .re = random.float(T) * @as(T, 10.0) - @as(T, 5.0),
                .im = random.float(T) * @as(T, 10.0) - @as(T, 5.0),
            };
            data[i] = original[i];
        }

        try fftInPlace(T, allocator, data);
        try ifftInPlace(T, allocator, data);

        for (0..size) |i| {
            try expectApproxEqRel(original[i].re, data[i].re, tol);
            try expectApproxEqRel(original[i].im, data[i].im, tol);
        }
    }
}

test "FFT validation: FFT-IFFT round trip f32" { try testRoundTripGeneric(f32); }
test "FFT validation: FFT-IFFT round trip f64" { try testRoundTripGeneric(f64); }

// 测试6: 验证线性性质
fn testLinearityGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tol = tolerance(T);

    const size: usize = 64;
    const a: T = @as(T, 2.5);
    const b: T = @as(T, -1.7);

    var prng = std.Random.DefaultPrng.init(99999);
    const random = prng.random();

    var x = try allocator.alloc(std.math.Complex(T), size);
    defer allocator.free(x);
    var y = try allocator.alloc(std.math.Complex(T), size);
    defer allocator.free(y);
    var combined = try allocator.alloc(std.math.Complex(T), size);
    defer allocator.free(combined);

    for (0..size) |i| {
        x[i] = std.math.Complex(T){ .re = random.float(T) * @as(T, 2.0) - @as(T, 1.0), .im = @as(T, 0.0) };
        y[i] = std.math.Complex(T){ .re = random.float(T) * @as(T, 2.0) - @as(T, 1.0), .im = @as(T, 0.0) };
        combined[i] = std.math.Complex(T){
            .re = a * x[i].re + b * y[i].re,
            .im = a * x[i].im + b * y[i].im,
        };
    }

    try fftInPlace(T, allocator, combined);
    try fftInPlace(T, allocator, x);
    try fftInPlace(T, allocator, y);

    for (0..size) |k| {
        const expected_re = a * x[k].re + b * y[k].re;
        const expected_im = a * x[k].im + b * y[k].im;
        try expectApproxEqRel(expected_re, combined[k].re, tol);
        try expectApproxEqRel(expected_im, combined[k].im, tol);
    }
}

test "FFT validation: Linearity property f32" { try testLinearityGeneric(f32); }
test "FFT validation: Linearity property f64" { try testLinearityGeneric(f64); }

// 测试7: 验证共轭对称性
fn testConjugateSymGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tol = tolerance(T);

    const sizes = [_]usize{ 8, 16, 64, 256 };
    var prng = std.Random.DefaultPrng.init(11111);
    const random = prng.random();

    for (sizes) |size| {
        var data = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(data);

        for (0..size) |i| {
            data[i] = std.math.Complex(T){ .re = random.float(T) * @as(T, 2.0) - @as(T, 1.0), .im = @as(T, 0.0) };
        }

        try fftInPlace(T, allocator, data);

        for (1..size / 2) |k| {
            const k_conj = size - k;
            try expectApproxEqRel(data[k].re, data[k_conj].re, tol);
            try expectApproxEqRel(data[k].im, -data[k_conj].im, tol);
        }
    }
}

test "FFT validation: Conjugate symmetry f32" { try testConjugateSymGeneric(f32); }
test "FFT validation: Conjugate symmetry f64" { try testConjugateSymGeneric(f64); }
