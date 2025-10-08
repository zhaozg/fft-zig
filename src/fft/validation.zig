//! Comprehensive FFT validation tests
//! 全面的FFT正确性验证测试

const std = @import("std");
const math = std.math;
const Complex = @import("types.zig").Complex;
const fftInPlace = @import("base.zig").fftInPlaceBase;
const ifft_module = @import("ifft.zig");
const ifftInPlace = ifft_module.ifftInPlace;

const expectApproxEqRel = std.testing.expectApproxEqRel;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;
const expect = std.testing.expect;
const TEST_TOLERANCE = 1e-10;

// 测试1: 验证单位冲激响应
// Unit impulse: FFT([1,0,0,...,0]) should be [1,1,1,...,1]
test "FFT validation: Unit impulse response" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const sizes = [_]usize{ 4, 8, 16, 64, 256, 1024 };
    for (sizes) |size| {
        var data = try allocator.alloc(Complex, size);
        defer allocator.free(data);

        // 创建单位冲激
        data[0] = Complex{ .re = 1.0, .im = 0.0 };
        for (1..size) |i| {
            data[i] = Complex{ .re = 0.0, .im = 0.0 };
        }

        try fftInPlace(allocator, data);

        // FFT([1,0,0,...,0]) = [1,1,1,...,1]
        for (data) |v| {
            try expectApproxEqRel(@as(f64, 1.0), v.re, TEST_TOLERANCE);
            try expectApproxEqAbs(@as(f64, 0.0), v.im, TEST_TOLERANCE);
        }
    }
}

// 测试2: 验证直流分量
// DC component: FFT([c,c,c,...,c]) should have all energy in bin 0
test "FFT validation: DC component" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const sizes = [_]usize{ 4, 8, 16, 64, 256 };
    for (sizes) |size| {
        var data = try allocator.alloc(Complex, size);
        defer allocator.free(data);

        const dc_value: f64 = 3.14159;
        for (0..size) |i| {
            data[i] = Complex{ .re = dc_value, .im = 0.0 };
        }

        try fftInPlace(allocator, data);

        // 直流分量应该在bin 0，值为 N * dc_value
        const expected_dc = @as(f64, @floatFromInt(size)) * dc_value;
        try expectApproxEqRel(expected_dc, data[0].re, TEST_TOLERANCE);
        try expectApproxEqAbs(@as(f64, 0.0), data[0].im, TEST_TOLERANCE);

        // 其他bins应该接近0
        for (1..size) |k| {
            const mag = @sqrt(data[k].re * data[k].re + data[k].im * data[k].im);
            try expect(mag < 1e-8);
        }
    }
}

// 测试3: 验证单频正弦波
// Single frequency: FFT(sin(2πk₀n/N)) should have peaks at k=k₀ and k=N-k₀
test "FFT validation: Single frequency sine wave" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const size: usize = 64;
    const freq: usize = 5; // 测试频率

    var data = try allocator.alloc(Complex, size);
    defer allocator.free(data);

    // 生成正弦波: sin(2π*freq*n/N)
    for (0..size) |n| {
        const angle = 2.0 * math.pi * @as(f64, @floatFromInt(freq * n)) / @as(f64, @floatFromInt(size));
        data[n] = Complex{ .re = @sin(angle), .im = 0.0 };
    }

    try fftInPlace(allocator, data);

    // 找到最大幅值的频率bin
    var max_mag: f64 = 0.0;
    var max_bin: usize = 0;
    for (0..size / 2) |k| {
        const mag = @sqrt(data[k].re * data[k].re + data[k].im * data[k].im);
        if (mag > max_mag) {
            max_mag = mag;
            max_bin = k;
        }
    }

    // 最大峰值应该在预期频率
    try expect(max_bin == freq);

    // 峰值幅度应该约为 N/2 (对于sin)
    const expected_peak = @as(f64, @floatFromInt(size)) / 2.0;
    try expectApproxEqRel(expected_peak, max_mag, 0.01);
}

// 测试4: 验证Parseval定理
// Parseval's theorem: Σ|x[n]|² = (1/N)·Σ|X[k]|²
test "FFT validation: Parseval's theorem" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const sizes = [_]usize{ 4, 8, 16, 64, 256, 1024 };
    var prng = std.Random.DefaultPrng.init(12345);
    const random = prng.random();

    for (sizes) |size| {
        var data = try allocator.alloc(Complex, size);
        defer allocator.free(data);

        // 生成随机信号
        var time_energy: f64 = 0.0;
        for (0..size) |i| {
            const val = random.float(f64) * 2.0 - 1.0;
            data[i] = Complex{ .re = val, .im = 0.0 };
            time_energy += val * val;
        }

        try fftInPlace(allocator, data);

        // 计算频域能量
        var freq_energy: f64 = 0.0;
        for (data) |v| {
            freq_energy += v.re * v.re + v.im * v.im;
        }

        // Parseval定理: time_energy = freq_energy / N
        const expected_freq_energy = time_energy * @as(f64, @floatFromInt(size));
        try expectApproxEqRel(expected_freq_energy, freq_energy, 0.01);
    }
}

// 测试5: 验证FFT-IFFT往返
// Round-trip: IFFT(FFT(x)) should equal x
test "FFT validation: FFT-IFFT round trip" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const sizes = [_]usize{ 4, 8, 16, 64, 256, 1024 };
    var prng = std.Random.DefaultPrng.init(54321);
    const random = prng.random();

    for (sizes) |size| {
        var original = try allocator.alloc(Complex, size);
        defer allocator.free(original);
        var data = try allocator.alloc(Complex, size);
        defer allocator.free(data);

        // 生成随机信号
        for (0..size) |i| {
            original[i] = Complex{
                .re = random.float(f64) * 10.0 - 5.0,
                .im = random.float(f64) * 10.0 - 5.0,
            };
            data[i] = original[i];
        }

        // FFT
        try fftInPlace(allocator, data);

        // IFFT
        try ifftInPlace(allocator, data);

        // 验证恢复到原始值
        for (0..size) |i| {
            try expectApproxEqRel(original[i].re, data[i].re, TEST_TOLERANCE);
            try expectApproxEqRel(original[i].im, data[i].im, TEST_TOLERANCE);
        }
    }
}

// 测试6: 验证线性性质
// Linearity: FFT(a·x + b·y) = a·FFT(x) + b·FFT(y)
test "FFT validation: Linearity property" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const size: usize = 64;
    const a: f64 = 2.5;
    const b: f64 = -1.7;

    var prng = std.Random.DefaultPrng.init(99999);
    const random = prng.random();

    var x = try allocator.alloc(Complex, size);
    defer allocator.free(x);
    var y = try allocator.alloc(Complex, size);
    defer allocator.free(y);
    var combined = try allocator.alloc(Complex, size);
    defer allocator.free(combined);

    // 生成随机信号
    for (0..size) |i| {
        x[i] = Complex{ .re = random.float(f64) * 2.0 - 1.0, .im = 0.0 };
        y[i] = Complex{ .re = random.float(f64) * 2.0 - 1.0, .im = 0.0 };
        combined[i] = Complex{
            .re = a * x[i].re + b * y[i].re,
            .im = a * x[i].im + b * y[i].im,
        };
    }

    // FFT(combined)
    try fftInPlace(allocator, combined);

    // FFT(x) and FFT(y)
    try fftInPlace(allocator, x);
    try fftInPlace(allocator, y);

    // 验证: FFT(a·x + b·y) = a·FFT(x) + b·FFT(y)
    for (0..size) |k| {
        const expected_re = a * x[k].re + b * y[k].re;
        const expected_im = a * x[k].im + b * y[k].im;
        try expectApproxEqRel(expected_re, combined[k].re, TEST_TOLERANCE);
        try expectApproxEqRel(expected_im, combined[k].im, TEST_TOLERANCE);
    }
}

// 测试7: 验证共轭对称性（实信号）
// Conjugate symmetry: For real input, X[k] = X*[N-k]
test "FFT validation: Conjugate symmetry for real signals" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const sizes = [_]usize{ 8, 16, 64, 256 };
    var prng = std.Random.DefaultPrng.init(11111);
    const random = prng.random();

    for (sizes) |size| {
        var data = try allocator.alloc(Complex, size);
        defer allocator.free(data);

        // 生成实信号
        for (0..size) |i| {
            data[i] = Complex{ .re = random.float(f64) * 2.0 - 1.0, .im = 0.0 };
        }

        try fftInPlace(allocator, data);

        // 验证共轭对称性: X[k] = conj(X[N-k])
        for (1..size / 2) |k| {
            const k_conj = size - k;
            try expectApproxEqRel(data[k].re, data[k_conj].re, TEST_TOLERANCE);
            try expectApproxEqRel(data[k].im, -data[k_conj].im, TEST_TOLERANCE);
        }
    }
}
