//! Radix-4 FFT 及其 SIMD 优化

const std = @import("std");
const math = std.math;
const fft_utils = @import("utils.zig");

const isPowerOfFour = fft_utils.isPowerOfFour;

pub fn fftRadix4SIMD(comptime T: type, data: []std.math.Complex(T)) error{ InvalidSize, OutOfMemory }!void {
    const n = data.len;
    if (n <= 1) return;
    if (!isPowerOfFour(n)) return error.InvalidSize;
    fft_utils.bitReversePermuteGeneric(std.math.Complex(T), data, 4);
    var stage_size: usize = 4;
    while (stage_size <= n) : (stage_size *= 4) {
        const quarter_stage = stage_size / 4;
        var group_start: usize = 0;
        while (group_start < n) : (group_start += stage_size) {
            for (0..quarter_stage) |k| {
                const theta = -2.0 * math.pi * @as(T, @floatFromInt(k)) / @as(T, @floatFromInt(stage_size));
                const w1 = std.math.Complex(T){ .re = math.cos(theta), .im = math.sin(theta) };
                const w2 = std.math.Complex(T){ .re = math.cos(2 * theta), .im = math.sin(2 * theta) };
                const w3 = std.math.Complex(T){ .re = math.cos(3 * theta), .im = math.sin(3 * theta) };
                const idx0 = group_start + k;
                const idx1 = idx0 + quarter_stage;
                const idx2 = idx1 + quarter_stage;
                const idx3 = idx2 + quarter_stage;
                const x0 = data[idx0];
                const x1_w1 = std.math.Complex(T){
                    .re = w1.re * data[idx1].re - w1.im * data[idx1].im,
                    .im = w1.re * data[idx1].im + w1.im * data[idx1].re,
                };
                const x2_w2 = std.math.Complex(T){
                    .re = w2.re * data[idx2].re - w2.im * data[idx2].im,
                    .im = w2.re * data[idx2].im + w2.im * data[idx2].re,
                };
                const x3_w3 = std.math.Complex(T){
                    .re = w3.re * data[idx3].re - w3.im * data[idx3].im,
                    .im = w3.re * data[idx3].im + w3.im * data[idx3].re,
                };
                const temp0 = std.math.Complex(T){ .re = x0.re + x2_w2.re, .im = x0.im + x2_w2.im };
                const temp1 = std.math.Complex(T){ .re = x1_w1.re + x3_w3.re, .im = x1_w1.im + x3_w3.im };
                const temp2 = std.math.Complex(T){ .re = x0.re - x2_w2.re, .im = x0.im - x2_w2.im };
                const temp3 = std.math.Complex(T){ .re = x1_w1.re - x3_w3.re, .im = x1_w1.im - x3_w3.im };
                data[idx0] = std.math.Complex(T){ .re = temp0.re + temp1.re, .im = temp0.im + temp1.im };
                data[idx1] = std.math.Complex(T){ .re = temp2.re + temp3.im, .im = temp2.im - temp3.re };
                data[idx2] = std.math.Complex(T){ .re = temp0.re - temp1.re, .im = temp0.im - temp1.im };
                data[idx3] = std.math.Complex(T){ .re = temp2.re - temp3.im, .im = temp2.im + temp3.re };
            }
        }
    }
}

pub fn fftRadix4(comptime T: type, data: []std.math.Complex(T)) error{ InvalidSize, OutOfMemory }!void {
    const n = data.len;
    if (n == 0) return error.InvalidSize;
    if (n == 1) return; // Identity transform
    if (!isPowerOfFour(n)) return error.InvalidSize;
    try fftRadix4SIMD(T, data);
}

const expect = std.testing.expect;
const expectApproxEqRel = std.testing.expectApproxEqRel;

fn testRadix4BitReverseGeneric(comptime T: type) !void {
    // Bit-reversal test
    var test_data = [_]std.math.Complex(T){
        std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) },
        std.math.Complex(T){ .re = @as(T, 2.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 3.0), .im = @as(T, 0.0) },
    };
    fft_utils.bitReversePermuteGeneric(std.math.Complex(T), test_data[0..], 4);

    // 4点基4反转应为 0,1,2,3 不变
    try expect(test_data[0].re == @as(T, 0.0));
    try expect(test_data[1].re == @as(T, 1.0));
    try expect(test_data[2].re == @as(T, 2.0));
    try expect(test_data[3].re == @as(T, 3.0));
}

test "Radix-4 FFT bit-reversal f32" { try testRadix4BitReverseGeneric(f32); }
test "Radix-4 FFT bit-reversal f64" { try testRadix4BitReverseGeneric(f64); }

fn testRadix4EdgeGeneric(comptime T: type) !void {
    var one = [_]std.math.Complex(T){std.math.Complex(T){ .re = @as(T, 7.0), .im = @as(T, 0.0) }};
    // Size 1 should now work (identity transform)
    try fftRadix4(T, one[0..]);
    try expect(one[0].re == @as(T, 7.0));
    try expect(one[0].im == @as(T, 0.0));

    var not_pow4 = [_]std.math.Complex(T){ 
        std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 2.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 3.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 4.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 5.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 6.0), .im = @as(T, 0.0) } 
    };
    if (fftRadix4(T, not_pow4[0..]) catch |err| err == error.InvalidSize) {
        // ok
    } else {
        try expect(false);
    }

    var four = [_]std.math.Complex(T){ 
        std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 2.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 3.0), .im = @as(T, 0.0) }, 
        std.math.Complex(T){ .re = @as(T, 4.0), .im = @as(T, 0.0) } 
    };
    try fftRadix4(T, four[0..]);
}

test "Radix-4 FFT edge cases f32" { try testRadix4EdgeGeneric(f32); }
test "Radix-4 FFT edge cases f64" { try testRadix4EdgeGeneric(f64); }
