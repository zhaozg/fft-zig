//! 混合基数与Bluestein算法及优化DFT

const std = @import("std");
const math = std.math;
const Complex = @import("types.zig").Complex;
const fft_radix2 = @import("fft_radix2.zig");
const fft_ifft = @import("ifft.zig");
const fft_utils = @import("utils.zig");
const fftRadix2 = fft_radix2.fftRadix2;
const ifftInPlace = fft_ifft.ifftInPlace;
const fftInPlaceBase = @import("base.zig").fftInPlaceBase;
// 循环依赖修复：fftInPlace 由主接口调度，不在此直接 import

pub fn fftMixedRadix(allocator: std.mem.Allocator, data: []Complex) error{ InvalidSize, OutOfMemory }!void {
    const n = data.len;
    if (n <= 1) return;
    if (isPowerOfTwo(n)) {
        return fftRadix2(data);
    }
    if (n < 1024)
        try optimizedDFTInPlace(data)
    else {
        try fftBluestein(allocator, data);
    }
}

pub fn fftBluestein(allocator: std.mem.Allocator, data: []Complex) !void {
    const n = data.len;
    if (n <= 1) return;
    const m = nextPowerOfTwo(2 * n - 1);
    var chirp = try allocator.alloc(Complex, n);
    defer allocator.free(chirp);
    for (0..n) |k| {
        const angle = math.pi * @as(f64, @floatFromInt(k * k)) / @as(f64, @floatFromInt(n));
        chirp[k] = Complex{
            .re = math.cos(angle),
            .im = -math.sin(angle),
        };
    }
    var a = try allocator.alloc(Complex, m);
    defer allocator.free(a);
    var b = try allocator.alloc(Complex, m);
    defer allocator.free(b);
    for (0..n) |k| {
        a[k] = data[k].mul(chirp[k]);
    }
    for (n..m) |k| {
        a[k] = Complex{ .re = 0.0, .im = 0.0 };
    }
    for (0..n) |k| {
        const angle = math.pi * @as(f64, @floatFromInt(k * k)) / @as(f64, @floatFromInt(n));
        b[k] = Complex{
            .re = math.cos(angle),
            .im = math.sin(angle),
        };
    }
    // Mirror chirp values for circular convolution
    for (1..n) |k| {
        b[m - k] = b[k];
    }
    for (n..(m - n + 1)) |k| {
        b[k] = Complex{ .re = 0.0, .im = 0.0 };
    }
    try fftInPlaceBase(allocator, a);
    try fftInPlaceBase(allocator, b);
    var c = try allocator.alloc(Complex, m);
    defer allocator.free(c);
    for (0..m) |k| {
        c[k] = a[k].mul(b[k]);
    }
    try ifftInPlace(allocator, c);
    for (0..n) |k| {
        data[k] = c[k].mul(chirp[k]);
    }
}

pub fn optimizedDFTInPlace(data: []Complex) error{ InvalidSize, OutOfMemory }!void {
    const n = data.len;
    if (n <= 1) return;
    const temp = try std.heap.page_allocator.alloc(Complex, n);
    defer std.heap.page_allocator.free(temp);
    const twiddle = try std.heap.page_allocator.alloc(Complex, n);
    defer std.heap.page_allocator.free(twiddle);
    for (0..n) |j| {
        const angle = -2.0 * math.pi * @as(f64, @floatFromInt(j)) / @as(f64, @floatFromInt(n));
        twiddle[j] = Complex{
            .re = math.cos(angle),
            .im = math.sin(angle),
        };
    }
    for (0..n) |k| {
        temp[k] = Complex{ .re = 0.0, .im = 0.0 };
        for (0..n) |j| {
            const twiddle_idx = (k * j) % n;
            const mult_result = Complex{
                .re = data[j].re * twiddle[twiddle_idx].re - data[j].im * twiddle[twiddle_idx].im,
                .im = data[j].re * twiddle[twiddle_idx].im + data[j].im * twiddle[twiddle_idx].re,
            };
            temp[k] = Complex{
                .re = temp[k].re + mult_result.re,
                .im = temp[k].im + mult_result.im,
            };
        }
    }
    @memcpy(data, temp);
}

fn isPowerOfTwo(n: usize) bool {
    return n > 0 and (n & (n - 1)) == 0;
}
fn nextPowerOfTwo(n: usize) usize {
    if (isPowerOfTwo(n)) return n;
    var power: usize = 1;
    while (power < n) power <<= 1;
    return power;
}

const expectApproxEqRel = std.testing.expectApproxEqRel;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;
const expect = std.testing.expect;
const TEST_TOLERANCE = 1e-10;

test "Mixed radix and DFT correctness" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // 8-point FFT vs DFT
    const test_data = [_]f64{ 1.0, 2.0, 1.0, -1.0, 1.5, 0.5, -0.5, 2.5 };
    var fft_data = try allocator.alloc(Complex, 8);
    defer allocator.free(fft_data);
    for (0..8) |i| {
        fft_data[i] = Complex{ .re = test_data[i], .im = 0.0 };
    }
    try fftMixedRadix(allocator, fft_data);

    var dft_data = try allocator.alloc(Complex, 8);
    defer allocator.free(dft_data);
    for (0..8) |k| {
        dft_data[k] = Complex{ .re = 0.0, .im = 0.0 };
        for (0..8) |n| {
            const angle = -2.0 * std.math.pi * @as(f64, @floatFromInt(k * n)) / 8.0;
            const w = Complex{ .re = @cos(angle), .im = @sin(angle) };
            const input_val = Complex{ .re = test_data[n], .im = 0.0 };
            dft_data[k] = Complex{
                .re = dft_data[k].re + input_val.re * w.re - input_val.im * w.im,
                .im = dft_data[k].im + input_val.re * w.im + input_val.im * w.re,
            };
        }
        // Forward DFT does NOT normalize
    }
    for (0..8) |i| {
        try expectApproxEqRel(fft_data[i].re, dft_data[i].re, TEST_TOLERANCE);
        try expectApproxEqAbs(fft_data[i].im, dft_data[i].im, TEST_TOLERANCE);
    }

    // Non-power-of-2 mixed radix
    {
        const size = 15; // 3 * 5
        var input = try allocator.alloc(Complex, size);
        defer allocator.free(input);

        for (0..size) |i| {
            input[i] = Complex{ .re = @as(f64, @floatFromInt(i)), .im = 0.0 };
        }

        try fftMixedRadix(allocator, input);

        const expected_dc = @as(f64, @floatFromInt((size - 1) * size / 2));
        try expectApproxEqRel(expected_dc, input[0].re, TEST_TOLERANCE);
    }
}

test "Mixed radix edge cases" {
    var empty: [0]Complex = undefined;
    try fftMixedRadix(std.heap.page_allocator, empty[0..]); // 应不报错

    var one = [_]Complex{Complex{ .re = 7.0, .im = 0.0 }};
    try fftMixedRadix(std.heap.page_allocator, one[0..]);
    try expectApproxEqRel(one[0].re, 7.0, TEST_TOLERANCE);

    var not_pow2 = [_]Complex{ Complex{ .re = 1.0, .im = 0.0 }, Complex{ .re = 2.0, .im = 0.0 }, Complex{ .re = 3.0, .im = 0.0 } };
    try fftMixedRadix(std.heap.page_allocator, not_pow2[0..]);
    // 可加断言检查输出
}
