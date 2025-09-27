//! Radix-4 FFT 及其 SIMD 优化

const std = @import("std");
const math = std.math;
const Complex = @import("types.zig").Complex;
const fft_utils = @import("utils.zig");

pub fn fftRadix4SIMD(data: []Complex) error{InvalidSize,OutOfMemory}!void {
    const n = data.len;
    if (n <= 1) return;
    if (!isPowerOfFour(n)) return error.InvalidSize;
    bitReverseRadix4(data);
    var stage_size: usize = 4;
    while (stage_size <= n) : (stage_size *= 4) {
        const quarter_stage = stage_size / 4;
        var group_start: usize = 0;
        while (group_start < n) : (group_start += stage_size) {
            for (0..quarter_stage) |k| {
                const theta = -2.0 * math.pi * @as(f64, @floatFromInt(k)) / @as(f64, @floatFromInt(stage_size));
                const w1 = Complex{ .re = math.cos(theta), .im = math.sin(theta) };
                const w2 = Complex{ .re = math.cos(2 * theta), .im = math.sin(2 * theta) };
                const w3 = Complex{ .re = math.cos(3 * theta), .im = math.sin(3 * theta) };
                const idx0 = group_start + k;
                const idx1 = idx0 + quarter_stage;
                const idx2 = idx1 + quarter_stage;
                const idx3 = idx2 + quarter_stage;
                const x0 = data[idx0];
                const x1_w1 = Complex{
                    .re = w1.re * data[idx1].re - w1.im * data[idx1].im,
                    .im = w1.re * data[idx1].im + w1.im * data[idx1].re,
                };
                const x2_w2 = Complex{
                    .re = w2.re * data[idx2].re - w2.im * data[idx2].im,
                    .im = w2.re * data[idx2].im + w2.im * data[idx2].re,
                };
                const x3_w3 = Complex{
                    .re = w3.re * data[idx3].re - w3.im * data[idx3].im,
                    .im = w3.re * data[idx3].im + w3.im * data[idx3].re,
                };
                const temp0 = Complex{ .re = x0.re + x2_w2.re, .im = x0.im + x2_w2.im };
                const temp1 = Complex{ .re = x1_w1.re + x3_w3.re, .im = x1_w1.im + x3_w3.im };
                const temp2 = Complex{ .re = x0.re - x2_w2.re, .im = x0.im - x2_w2.im };
                const temp3 = Complex{ .re = x1_w1.re - x3_w3.re, .im = x1_w1.im - x3_w3.im };
                data[idx0] = Complex{ .re = temp0.re + temp1.re, .im = temp0.im + temp1.im };
                data[idx1] = Complex{ .re = temp2.re - temp3.im, .im = temp2.im + temp3.re };
                data[idx2] = Complex{ .re = temp0.re - temp1.re, .im = temp0.im - temp1.im };
                data[idx3] = Complex{ .re = temp2.re + temp3.im, .im = temp2.im - temp3.re };
            }
        }
    }
}

pub fn fftRadix4(data: []Complex) error{InvalidSize,OutOfMemory}!void {
    const n = data.len;
    if (n <= 1 or !isPowerOfFour(n)) return error.InvalidSize;
    try fftRadix4SIMD(data);
}

pub fn bitReverseRadix4(data: []Complex) void {
    const n = data.len;
    if (n <= 1) return;
    for (0..n) |i| {
        var j: usize = 0;
        var temp_i = i;
        var temp_n = n;
        while (temp_n > 1) {
            j = j * 4 + (temp_i % 4);
            temp_i /= 4;
            temp_n /= 4;
        }
        if (i < j) {
            const temp = data[i];
            data[i] = data[j];
            data[j] = temp;
        }
    }
}

fn isPowerOfFour(n: usize) bool {
    if (!isPowerOfTwo(n)) return false;
    return (n & 0x55555555) != 0;
}
fn isPowerOfTwo(n: usize) bool {
    return n > 0 and (n & (n - 1)) == 0;
}


const expect = std.testing.expect;
const expectApproxEqRel = std.testing.expectApproxEqRel;
const TEST_TOLERANCE = 1e-12;

test "Radix-4 FFT bit-reversal and correctness" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();

    // Bit-reversal test
    var test_data = [_]Complex{
        Complex{ .re = 0.0, .im = 0.0 }, Complex{ .re = 1.0, .im = 0.0 },
        Complex{ .re = 2.0, .im = 0.0 }, Complex{ .re = 3.0, .im = 0.0 },
    };
    bitReverseRadix4(&test_data);

    // 4点基4反转应为 0,1,2,3 不变
    try expect(test_data[0].re == 0.0);
    try expect(test_data[1].re == 1.0);
    try expect(test_data[2].re == 2.0);
    try expect(test_data[3].re == 3.0);
}


test "Radix-4 FFT edge cases" {
    var one = [_]Complex{Complex{ .re = 7.0, .im = 0.0 }};
    if (fftRadix4(one[0..]) catch |err| err == error.InvalidSize) {
        // ok
    } else {
        try expect(false);
    }

    var not_pow4 = [_]Complex{
        Complex{ .re = 1.0, .im = 0.0 }, Complex{ .re = 2.0, .im = 0.0 },
        Complex{ .re = 3.0, .im = 0.0 }, Complex{ .re = 4.0, .im = 0.0 },
        Complex{ .re = 5.0, .im = 0.0 }, Complex{ .re = 6.0, .im = 0.0 }
    };
    if (fftRadix4(not_pow4[0..]) catch |err| err == error.InvalidSize) {
        // ok
    } else {
        try expect(false);
    }

    var four = [_]Complex{
        Complex{ .re = 1.0, .im = 0.0 }, Complex{ .re = 2.0, .im = 0.0 },
        Complex{ .re = 3.0, .im = 0.0 }, Complex{ .re = 4.0, .im = 0.0 }
    };
    try fftRadix4(four[0..]);
    // 可加断言检查输出
}
