//! 混合基数与Bluestein算法及优化DFT

const std = @import("std");
const math = std.math;
const fft_radix2 = @import("fft_radix2.zig");
const fft_ifft = @import("ifft.zig");
// 循环依赖修复：fftInPlace 由主接口调度，不在此直接 import

pub fn fftMixedRadix(comptime T: type, allocator: std.mem.Allocator, data: []std.math.Complex(T)) error{ InvalidSize, OutOfMemory }!void {
    const n = data.len;
    if (n <= 1) return;
    if (isPowerOfTwo(n)) {
        return fft_radix2.fftRadix2(T, data);
    }
    if (n < 1024)
        try optimizedDFTInPlace(T, data)
    else {
        try fftBluestein(T, allocator, data);
    }
}

pub fn fftBluestein(comptime T: type, allocator: std.mem.Allocator, data: []std.math.Complex(T)) !void {
    const n = data.len;
    if (n <= 1) return;
    const m = nextPowerOfTwo(2 * n - 1);
    var chirp = try allocator.alloc(std.math.Complex(T), n);
    defer allocator.free(chirp);
    for (0..n) |k| {
        const angle = math.pi * @as(T, @floatFromInt(k * k)) / @as(T, @floatFromInt(n));
        chirp[k] = std.math.Complex(T){
            .re = math.cos(angle),
            .im = -math.sin(angle),
        };
    }
    var a = try allocator.alloc(std.math.Complex(T), m);
    defer allocator.free(a);
    var b = try allocator.alloc(std.math.Complex(T), m);
    defer allocator.free(b);
    for (0..n) |k| {
        a[k] = data[k].mul(chirp[k]);
    }
    for (n..m) |k| {
        a[k] = std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) };
    }
    for (0..n) |k| {
        const angle = math.pi * @as(T, @floatFromInt(k * k)) / @as(T, @floatFromInt(n));
        b[k] = std.math.Complex(T){
            .re = math.cos(angle),
            .im = math.sin(angle),
        };
    }
    // Mirror chirp values for circular convolution
    for (1..n) |k| {
        b[m - k] = b[k];
    }
    for (n..(m - n + 1)) |k| {
        b[k] = std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) };
    }
    try @import("base.zig").fftInPlaceBase(T, allocator, a);
    try @import("base.zig").fftInPlaceBase(T, allocator, b);
    var c = try allocator.alloc(std.math.Complex(T), m);
    defer allocator.free(c);
    for (0..m) |k| {
        c[k] = a[k].mul(b[k]);
    }
    try fft_ifft.ifftInPlace(T, allocator, c);
    for (0..n) |k| {
        data[k] = c[k].mul(chirp[k]);
    }
}

pub fn optimizedDFTInPlace(comptime T: type, data: []std.math.Complex(T)) error{ InvalidSize, OutOfMemory }!void {
    const n = data.len;
    if (n <= 1) return;
    const temp = try std.heap.page_allocator.alloc(std.math.Complex(T), n);
    defer std.heap.page_allocator.free(temp);
    const twiddle = try std.heap.page_allocator.alloc(std.math.Complex(T), n);
    defer std.heap.page_allocator.free(twiddle);
    for (0..n) |j| {
        const angle = -2.0 * math.pi * @as(T, @floatFromInt(j)) / @as(T, @floatFromInt(n));
        twiddle[j] = std.math.Complex(T){
            .re = math.cos(angle),
            .im = math.sin(angle),
        };
    }
    for (0..n) |k| {
        temp[k] = std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) };
        for (0..n) |j| {
            const twiddle_idx = (k * j) % n;
            const mult_result = std.math.Complex(T){
                .re = data[j].re * twiddle[twiddle_idx].re - data[j].im * twiddle[twiddle_idx].im,
                .im = data[j].re * twiddle[twiddle_idx].im + data[j].im * twiddle[twiddle_idx].re,
            };
            temp[k] = std.math.Complex(T){
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

fn testMixedRadixGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tolerance: T = if (T == f32) @as(T, 1e-5) else @as(T, 1e-10);

    // 8-point FFT vs DFT
    const test_data = [_]T{ @as(T, 1.0), @as(T, 2.0), @as(T, 1.0), @as(T, -1.0), @as(T, 1.5), @as(T, 0.5), @as(T, -0.5), @as(T, 2.5) };
    var fft_data = try allocator.alloc(std.math.Complex(T), 8);
    defer allocator.free(fft_data);
    for (0..8) |i| {
        fft_data[i] = std.math.Complex(T){ .re = test_data[i], .im = @as(T, 0.0) };
    }
    try fftMixedRadix(T, allocator, fft_data);

    var dft_data = try allocator.alloc(std.math.Complex(T), 8);
    defer allocator.free(dft_data);
    for (0..8) |k| {
        dft_data[k] = std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) };
        for (0..8) |n| {
            const angle = -2.0 * std.math.pi * @as(T, @floatFromInt(k * n)) / @as(T, 8.0);
            const w = std.math.Complex(T){ .re = @cos(angle), .im = @sin(angle) };
            const input_val = std.math.Complex(T){ .re = test_data[n], .im = @as(T, 0.0) };
            dft_data[k] = std.math.Complex(T){
                .re = dft_data[k].re + input_val.re * w.re - input_val.im * w.im,
                .im = dft_data[k].im + input_val.re * w.im + input_val.im * w.re,
            };
        }
        // Forward DFT does NOT normalize
    }
    for (0..8) |i| {
        try expectApproxEqRel(fft_data[i].re, dft_data[i].re, tolerance);
        try expectApproxEqAbs(fft_data[i].im, dft_data[i].im, tolerance);
    }

    // Non-power-of-2 mixed radix
    {
        const size = 15; // 3 * 5
        var input = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(input);

        for (0..size) |i| {
            input[i] = std.math.Complex(T){ .re = @as(T, @floatFromInt(i)), .im = @as(T, 0.0) };
        }

        try fftMixedRadix(T, allocator, input);

        const expected_dc = @as(T, @floatFromInt((size - 1) * size / 2));
        try expectApproxEqRel(expected_dc, input[0].re, tolerance);
    }
}

test "Mixed radix and DFT correctness f32" { try testMixedRadixGeneric(f32); }
test "Mixed radix and DFT correctness f64" { try testMixedRadixGeneric(f64); }

fn testMixedRadixEdgeGeneric(comptime T: type) !void {
    const tolerance: T = if (T == f32) @as(T, 1e-6) else @as(T, 1e-12);
    var empty: [0]std.math.Complex(T) = undefined;
    try fftMixedRadix(T, std.heap.page_allocator, empty[0..]); // 应不报错

    var one = [_]std.math.Complex(T){std.math.Complex(T){ .re = @as(T, 7.0), .im = @as(T, 0.0) }};
    try fftMixedRadix(T, std.heap.page_allocator, one[0..]);
    try expectApproxEqRel(one[0].re, @as(T, 7.0), tolerance);

    var not_pow2 = [_]std.math.Complex(T){ 
        std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 2.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 3.0), .im = @as(T, 0.0) } 
    };
    try fftMixedRadix(T, std.heap.page_allocator, not_pow2[0..]);
}

test "Mixed radix edge cases f32" { try testMixedRadixEdgeGeneric(f32); }
test "Mixed radix edge cases f64" { try testMixedRadixEdgeGeneric(f64); }
