//! Radix-2 FFT 及其 SIMD 优化

const std = @import("std");
const math = std.math;
const Complex = @import("types.zig").Complex;
const VectorF64 = @import("types.zig").VectorF64;
const fft_utils = @import("utils.zig");

pub fn fftRadix2(data: []Complex) error{InvalidSize,OutOfMemory}!void {
    const n = data.len;
    if (n <= 1) return;
    if (!isPowerOfTwo(n)) return error.InvalidSize;
    bitReversePermute(data);
    var stage_size: usize = 2;
    while (stage_size <= n) : (stage_size *= 2) {
        const half_stage = stage_size / 2;
        const theta = -2.0 * math.pi / @as(f64, @floatFromInt(stage_size));
        var group_start: usize = 0;
        while (group_start < n) : (group_start += stage_size) {
            for (0..half_stage) |k| {
                const even_idx = group_start + k;
                const odd_idx = even_idx + half_stage;
                const angle = theta * @as(f64, @floatFromInt(k));
                const w = Complex{
                    .re = math.cos(angle),
                    .im = math.sin(angle),
                };
                const temp_re = w.re * data[odd_idx].re - w.im * data[odd_idx].im;
                const temp_im = w.re * data[odd_idx].im + w.im * data[odd_idx].re;
                data[odd_idx].re = data[even_idx].re - temp_re;
                data[odd_idx].im = data[even_idx].im - temp_im;
                data[even_idx].re = data[even_idx].re + temp_re;
                data[even_idx].im = data[even_idx].im + temp_im;
            }
        }
    }
}

pub fn fftRadix2SIMD(data: []Complex) error{InvalidSize,OutOfMemory}!void {
    const n = data.len;
    if (n <= 1) return;
    if (!isPowerOfTwo(n)) return error.InvalidSize;
    bitReversePermuteSIMD(data);
    var stage_size: usize = 2;
    while (stage_size <= n) : (stage_size *= 2) {
        const half_stage = stage_size / 2;
        const theta = -2.0 * math.pi / @as(f64, @floatFromInt(stage_size));
        var group_start: usize = 0;
        while (group_start < n) : (group_start += stage_size) {
            var k: usize = 0;
            while (k + 3 < half_stage) : (k += 4) {
                const k_vec = VectorF64{ @as(f64, @floatFromInt(k)), @as(f64, @floatFromInt(k + 1)), @as(f64, @floatFromInt(k + 2)), @as(f64, @floatFromInt(k + 3)) };
                const angles = k_vec * @as(VectorF64, @splat(theta));
                const cos_vals = @cos(angles);
                const sin_vals = @sin(angles);
                for (0..4) |i| {
                    const even_idx = group_start + k + i;
                    const odd_idx = even_idx + half_stage;
                    if (odd_idx >= n) break;
                    const w_re = cos_vals[i];
                    const w_im = sin_vals[i];
                    const temp_re = w_re * data[odd_idx].re - w_im * data[odd_idx].im;
                    const temp_im = w_re * data[odd_idx].im + w_im * data[odd_idx].re;
                    data[odd_idx].re = data[even_idx].re - temp_re;
                    data[odd_idx].im = data[even_idx].im - temp_im;
                    data[even_idx].re = data[even_idx].re + temp_re;
                    data[even_idx].im = data[even_idx].im + temp_im;
                }
            }
            while (k < half_stage) : (k += 1) {
                const even_idx = group_start + k;
                const odd_idx = even_idx + half_stage;
                const angle = theta * @as(f64, @floatFromInt(k));
                const w = Complex{
                    .re = math.cos(angle),
                    .im = math.sin(angle),
                };
                const temp_re = w.re * data[odd_idx].re - w.im * data[odd_idx].im;
                const temp_im = w.re * data[odd_idx].im + w.im * data[odd_idx].re;
                data[odd_idx].re = data[even_idx].re - temp_re;
                data[odd_idx].im = data[even_idx].im - temp_im;
                data[even_idx].re = data[even_idx].re + temp_re;
                data[even_idx].im = data[even_idx].im + temp_im;
            }
        }
    }
}

pub fn bitReversePermuteSIMD(data: []Complex) void {
    const n = data.len;
    if (n <= 1) return;
    var j: usize = 0;
    for (1..n) |i| {
        var bit = n >> 1;
        while (j & bit != 0) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) {
            const temp = data[i];
            data[i] = data[j];
            data[j] = temp;
        }
    }
}

pub fn bitReversePermute(data: []Complex) void {
    bitReversePermuteSIMD(data);
}

fn isPowerOfTwo(n: usize) bool {
    return n > 0 and (n & (n - 1)) == 0;
}

const expectApproxEqRel = std.testing.expectApproxEqRel;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;
const expect = std.testing.expect;
const TEST_TOLERANCE = 1e-12;

test "Radix-2 FFT basic functionality" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    // Size 1 FFT
    {
        var data = [_]Complex{Complex{ .re = 42.0, .im = 0.0 }};
        try fftRadix2(data[0..]);
        try expectApproxEqRel(@as(f64, 42.0), data[0].re, TEST_TOLERANCE);
    }
    // Size 2 FFT
    {
        var data = [_]Complex{
            Complex{ .re = 1.0, .im = 0.0 },
            Complex{ .re = -1.0, .im = 0.0 },
        };
        try fftRadix2(data[0..]);
        try expectApproxEqRel(@as(f64, 0.0), data[0].re, TEST_TOLERANCE);
        try expectApproxEqAbs(@as(f64, 2.0), data[1].re, TEST_TOLERANCE);
    }
}

test "Radix-2 FFT edge cases" {
    var empty: [0]Complex = undefined;
    try fftRadix2(empty[0..]); // 应不报错

    var one = [_]Complex{Complex{ .re = 7.0, .im = 0.0 }};
    try fftRadix2(one[0..]);
    try expectApproxEqRel(one[0].re, 7.0, TEST_TOLERANCE);

    var not_pow2 = [_]Complex{Complex{ .re = 1.0, .im = 0.0 }, Complex{ .re = 2.0, .im = 0.0 }, Complex{ .re = 3.0, .im = 0.0 }};
    const result = fftRadix2(not_pow2[0..]);
    try expect(result == error.InvalidSize);
}
