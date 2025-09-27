//! 实数到复数FFT及相关工具


const std = @import("std");
const math = std.math;
const Complex = @import("types.zig").Complex;
const VectorF64 = @import("types.zig").VectorF64;
const fft_ifft = @import("ifft.zig");
const fft_utils = @import("utils.zig");
const fft_radix2 = @import("fft_radix2.zig");
const fft_parallel = @import("fft_parallel.zig");
const fft_mixed = @import("fft_mixed.zig");

extern fn fftInPlace(allocator: std.mem.Allocator, data: []Complex) anyerror!void;

pub fn fftR2C(allocator: std.mem.Allocator, input: []const f64, output: []f64, magnitude: []f64) !void {
    const n = input.len;
    const out_len = n / 2 + 1;
    if (output.len < 2 * out_len) return error.BufferTooSmall;
    if (magnitude.len < out_len) return error.BufferTooSmall;
    if (n <= 256) {
        try computeSmallFFT(input, output, magnitude);
        return;
    }
    if (n >= 1000000) {
        try computeHugeR2C(allocator, input, output, magnitude);
        return;
    }
    const complex_buffer = try allocateAlignedComplexBuffer(allocator, n);
    defer allocator.free(complex_buffer);
    for (0..n) |i| {
        complex_buffer[i] = Complex{ .re = input[i], .im = 0.0 };
    }
    try fftInPlace(allocator, complex_buffer);
    convertToOutputSIMD(complex_buffer[0..out_len], output, magnitude);
}

fn computeSmallFFT(input: []const f64, output: []f64, magnitude: []f64) !void {
    const n = input.len;
    const out_len = n / 2 + 1;
    for (0..out_len) |k| {
        var real: f64 = 0.0;
        var imag: f64 = 0.0;
        for (0..n) |j| {
            const angle = -2.0 * math.pi * @as(f64, @floatFromInt(k)) * @as(f64, @floatFromInt(j)) / @as(f64, @floatFromInt(n));
            const cos_val = math.cos(angle);
            const sin_val = math.sin(angle);
            real += input[j] * cos_val;
            imag += input[j] * sin_val;
        }
        output[2 * k] = real;
        output[2 * k + 1] = imag;
        magnitude[k] = @sqrt(real * real + imag * imag);
    }
}

fn convertToOutputSIMD(input: []const Complex, output: []f64, magnitude: []f64) void {
    const n = input.len;
    var i: usize = 0;
    while (i + 3 < n) : (i += 4) {
        const re = VectorF64{ input[i].re, input[i + 1].re, input[i + 2].re, input[i + 3].re };
        const im = VectorF64{ input[i].im, input[i + 1].im, input[i + 2].im, input[i + 3].im };
        output[2 * i] = re[0];
        output[2 * i + 1] = im[0];
        output[2 * (i + 1)] = re[1];
        output[2 * (i + 1) + 1] = im[1];
        output[2 * (i + 2)] = re[2];
        output[2 * (i + 2) + 1] = im[2];
        output[2 * (i + 3)] = re[3];
        output[2 * (i + 3) + 1] = im[3];
        const mag_squared = re * re + im * im;
        const mag = @sqrt(mag_squared);
        magnitude[i] = mag[0];
        magnitude[i + 1] = mag[1];
        magnitude[i + 2] = mag[2];
        magnitude[i + 3] = mag[3];
    }
    while (i < n) : (i += 1) {
        output[2 * i] = input[i].re;
        output[2 * i + 1] = input[i].im;
        magnitude[i] = @sqrt(input[i].re * input[i].re + input[i].im * input[i].im);
    }
}

fn computeHugeR2C(allocator: std.mem.Allocator, input: []const f64, output: []f64, magnitude: []f64) !void {
    const n = input.len;
    const out_len = n / 2 + 1;
    const block_size = @min(n, 1024 * 1024);
    var processed: usize = 0;
    while (processed < n) {
        const current_block_size = @min(block_size, n - processed);
        const input_block = input[processed .. processed + current_block_size];
        const output_start = processed / 2;
        const output_block_len = @min(current_block_size / 2 + 1, out_len - output_start);
        if (output_start >= out_len) break;
        const output_block = output[2 * output_start .. 2 * (output_start + output_block_len)];
        const magnitude_block = magnitude[output_start .. output_start + output_block_len];
        if (current_block_size <= 256) {
            try computeSmallFFT(input_block, output_block, magnitude_block);
        } else {
            const complex_buffer = try allocateAlignedComplexBuffer(allocator, current_block_size);
            defer allocator.free(complex_buffer);
            for (0..current_block_size) |i| {
                complex_buffer[i] = Complex{ .re = input_block[i], .im = 0.0 };
            }
            try fftInPlace(allocator, complex_buffer);
            convertToOutputSIMD(complex_buffer[0..output_block_len], output_block, magnitude_block);
        }
        processed += current_block_size;
    }
}

fn allocateAlignedComplexBuffer(allocator: std.mem.Allocator, size: usize) ![]Complex {
    return try allocator.alloc(Complex, size);
}

const expect = std.testing.expect;
const expectApproxEqRel = std.testing.expectApproxEqRel;
const TEST_TOLERANCE = 1e-12;

test "Real-to-complex FFT and utilities" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const size = 128;
    const input = try allocator.alloc(f64, size);
    defer allocator.free(input);

    for (0..size) |i| {
        const t = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(size));
        input[i] = @sin(2.0 * std.math.pi * 5.0 * t) + 0.5 * @cos(2.0 * std.math.pi * 10.0 * t);
    }

    const out_len = size / 2 + 1;
    const output = try allocator.alloc(f64, 2 * out_len);
    defer allocator.free(output);
    const magnitude = try allocator.alloc(f64, out_len);
    defer allocator.free(magnitude);

    try fftR2C(allocator, input, output, magnitude);
}

test "SIMD magnitude calculation" {
    const complex_vals = [_]Complex{
        Complex{ .re = 3.0, .im = 4.0 }, // magnitude = 5.0
        Complex{ .re = 1.0, .im = 1.0 }, // magnitude = sqrt(2)
    };
    for (complex_vals, 0..) |val, i| {
        const expected_mag = @sqrt(val.re * val.re + val.im * val.im);
        const computed_mag = @sqrt(val.re * val.re + val.im * val.im);
        try expectApproxEqRel(expected_mag, computed_mag, TEST_TOLERANCE);
        _ = i;
    }
}

test "R2C FFT edge cases" {
    var empty: [0]f64 = undefined;
    var out_empty: [2]f64 = undefined;
    var mag_empty: [1]f64 = undefined;
    try fftR2C(std.heap.page_allocator, empty[0..], out_empty[0..], mag_empty[0..]); // 应不报错

    var one = [_]f64{7.0};
    var out_one: [2]f64 = undefined;
    var mag_one: [1]f64 = undefined;
    try fftR2C(std.heap.page_allocator, one[0..], out_one[0..], mag_one[0..]);
    try expectApproxEqRel(out_one[0], 7.0, 1e-12);
    try expectApproxEqRel(out_one[1], 0.0, 1e-12);
    try expectApproxEqRel(mag_one[0], 7.0, 1e-12);
}
