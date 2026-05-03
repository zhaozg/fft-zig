//! Radix-2 FFT 及其 SIMD 优化

const std = @import("std");
const math = std.math;
const fft_utils = @import("utils.zig");
const isPowerOfTwo = fft_utils.isPowerOfTwo;

pub fn fftRadix2(comptime T: type, data: []std.math.Complex(T)) error{ InvalidSize, OutOfMemory }!void {
    const n = data.len;
    if (n <= 1) return;
    if (!isPowerOfTwo(n)) return error.InvalidSize;
    fft_utils.bitReversePermuteGeneric(std.math.Complex(T), data, 2);
    var stage_size: usize = 2;
    while (stage_size <= n) : (stage_size *= 2) {
        const half_stage = stage_size / 2;
        // 递推法生成twiddle factor
        const theta = -2.0 * math.pi / @as(T, @floatFromInt(stage_size));
        const w_unit = std.math.Complex(T){ .re = math.cos(theta), .im = math.sin(theta) };
        var w = std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) };
        var group_start: usize = 0;
        while (group_start < n) : (group_start += stage_size) {
            w = std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) };
            for (0..half_stage) |k| {
                const even_idx = group_start + k;
                const odd_idx = even_idx + half_stage;
                const temp_re = w.re * data[odd_idx].re - w.im * data[odd_idx].im;
                const temp_im = w.re * data[odd_idx].im + w.im * data[odd_idx].re;
                data[odd_idx].re = data[even_idx].re - temp_re;
                data[odd_idx].im = data[even_idx].im - temp_im;
                data[even_idx].re = data[even_idx].re + temp_re;
                data[even_idx].im = data[even_idx].im + temp_im;
                // 递推更新w
                const next_re = w.re * w_unit.re - w.im * w_unit.im;
                const next_im = w.re * w_unit.im + w.im * w_unit.re;
                w.re = next_re;
                w.im = next_im;
            }
        }
    }
    // Forward FFT does NOT normalize - only IFFT should normalize
}

pub fn fftRadix2SIMD(comptime T: type, data: []std.math.Complex(T)) error{ InvalidSize, OutOfMemory }!void {
    const n = data.len;
    if (n <= 1) return;
    if (!isPowerOfTwo(n)) return error.InvalidSize;
    fft_utils.bitReversePermuteGeneric(std.math.Complex(T), data, 2);
    var stage_size: usize = 2;
    while (stage_size <= n) : (stage_size *= 2) {
        const half_stage = stage_size / 2;
        const theta = -2.0 * math.pi / @as(T, @floatFromInt(stage_size));
        const w_unit = std.math.Complex(T){ .re = math.cos(theta), .im = math.sin(theta) };
        var group_start: usize = 0;
        while (group_start < n) : (group_start += stage_size) {
            var w = std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) };
            var k: usize = 0;
            while (k + 3 < half_stage) : (k += 4) {
                // 递推生成4个twiddle
                var ws: [4]std.math.Complex(T) = undefined;
                ws[0] = w;
                for (1..4) |i| {
                    ws[i] = std.math.Complex(T){
                        .re = ws[i - 1].re * w_unit.re - ws[i - 1].im * w_unit.im,
                        .im = ws[i - 1].re * w_unit.im + ws[i - 1].im * w_unit.re,
                    };
                }
                for (0..4) |i| {
                    const even_idx = group_start + k + i;
                    const odd_idx = even_idx + half_stage;
                    if (odd_idx >= n) break;
                    const w_re = ws[i].re;
                    const w_im = ws[i].im;
                    const temp_re = w_re * data[odd_idx].re - w_im * data[odd_idx].im;
                    const temp_im = w_re * data[odd_idx].im + w_im * data[odd_idx].re;
                    data[odd_idx].re = data[even_idx].re - temp_re;
                    data[odd_idx].im = data[even_idx].im - temp_im;
                    data[even_idx].re = data[even_idx].re + temp_re;
                    data[even_idx].im = data[even_idx].im + temp_im;
                }
                // 递推到下一个w
                w = ws[3];
                const next_re = w.re * w_unit.re - w.im * w_unit.im;
                const next_im = w.re * w_unit.im + w.im * w_unit.re;
                w.re = next_re;
                w.im = next_im;
            }
            while (k < half_stage) : (k += 1) {
                const even_idx = group_start + k;
                const odd_idx = even_idx + half_stage;
                const temp_re = w.re * data[odd_idx].re - w.im * data[odd_idx].im;
                const temp_im = w.re * data[odd_idx].im + w.im * data[odd_idx].re;
                data[odd_idx].re = data[even_idx].re - temp_re;
                data[odd_idx].im = data[even_idx].im - temp_im;
                data[even_idx].re = data[even_idx].re + temp_re;
                data[even_idx].im = data[even_idx].im + temp_im;
                // 递推更新w
                const next_re = w.re * w_unit.re - w.im * w_unit.im;
                const next_im = w.re * w_unit.im + w.im * w_unit.re;
                w.re = next_re;
                w.im = next_im;
            }
        }
    }
    // Forward FFT does NOT normalize - only IFFT should normalize
}

const expectApproxEqRel = std.testing.expectApproxEqRel;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;
const expect = std.testing.expect;

fn testRadix2BasicGeneric(comptime T: type) !void {
    const tolerance: T = if (T == f32) @as(T, 1e-6) else @as(T, 1e-12);
    // Size 1 FFT
    {
        var data = [_]std.math.Complex(T){std.math.Complex(T){ .re = @as(T, 42.0), .im = @as(T, 0.0) }};
        try fftRadix2(T, data[0..]);
        try expectApproxEqRel(@as(T, 42.0), data[0].re, tolerance);
    }
    // Size 2 FFT - without normalization, [1, -1] -> [0, 2]
    {
        var data = [_]std.math.Complex(T){
            std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) },
            std.math.Complex(T){ .re = @as(T, -1.0), .im = @as(T, 0.0) },
        };
        try fftRadix2(T, data[0..]);
        try expectApproxEqRel(@as(T, 0.0), data[0].re, tolerance);
        try expectApproxEqAbs(@as(T, 2.0), data[1].re, tolerance);
    }
}

test "Radix-2 FFT basic functionality f32" { try testRadix2BasicGeneric(f32); }
test "Radix-2 FFT basic functionality f64" { try testRadix2BasicGeneric(f64); }

fn testRadix2EdgeGeneric(comptime T: type) !void {
    const tolerance: T = if (T == f32) @as(T, 1e-6) else @as(T, 1e-12);
    var empty: [0]std.math.Complex(T) = undefined;
    try fftRadix2(T, empty[0..]); // 应不报错

    var one = [_]std.math.Complex(T){std.math.Complex(T){ .re = @as(T, 7.0), .im = @as(T, 0.0) }};
    try fftRadix2(T, one[0..]);
    try expectApproxEqRel(one[0].re, @as(T, 7.0), tolerance);

    var not_pow2 = [_]std.math.Complex(T){ 
        std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 2.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 3.0), .im = @as(T, 0.0) } 
    };
    const result = fftRadix2(T, not_pow2[0..]);
    try expect(result == error.InvalidSize);
}

test "Radix-2 FFT edge cases f32" { try testRadix2EdgeGeneric(f32); }
test "Radix-2 FFT edge cases f64" { try testRadix2EdgeGeneric(f64); }
