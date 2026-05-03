//! 实数到复数FFT及相关工具

const std = @import("std");
const math = std.math;
const fft_types = @import("types.zig");

const fftInPlace = @import("fft_inplace.zig").fftInPlace;

pub fn fftR2C(comptime T: type, allocator: std.mem.Allocator, input: []const T, output: []T, magnitude: []T) !void {
    const n = input.len;
    const out_len = n / 2 + 1;
    if (output.len < 2 * out_len) return error.BufferTooSmall;
    if (magnitude.len < out_len) return error.BufferTooSmall;
    if (n <= 256) {
        try computeSmallFFT(T, input, output, magnitude);
        return;
    }
    if (n >= 1000000) {
        try computeHugeR2C(T, allocator, input, output, magnitude);
        return;
    }
    const complex_buffer = try allocator.alloc(std.math.Complex(T), n);
    defer allocator.free(complex_buffer);
    for (0..n) |i| {
        complex_buffer[i] = std.math.Complex(T){ .re = input[i], .im = @as(T, 0.0) };
    }
    try fftInPlace(T, allocator, complex_buffer);
    // Convert complex buffer to output format
    convertToOutput(T, complex_buffer[0..out_len], output, magnitude);
}

fn computeSmallFFT(comptime T: type, input: []const T, output: []T, magnitude: []T) !void {
    const n = input.len;
    const out_len = n / 2 + 1;
    for (0..out_len) |k| {
        var real: T = @as(T, 0.0);
        var imag: T = @as(T, 0.0);
        for (0..n) |j| {
            const angle = -2.0 * math.pi * @as(T, @floatFromInt(k)) * @as(T, @floatFromInt(j)) / @as(T, @floatFromInt(n));
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

fn convertToOutput(comptime T: type, input: []const std.math.Complex(T), output: []T, magnitude: []T) void {
    const n = input.len;
    for (0..n) |i| {
        output[2 * i] = input[i].re;
        output[2 * i + 1] = input[i].im;
        magnitude[i] = @sqrt(input[i].re * input[i].re + input[i].im * input[i].im);
    }
}

fn computeHugeR2C(comptime T: type, allocator: std.mem.Allocator, input: []const T, output: []T, magnitude: []T) !void {
    const n = input.len;
    const out_len = n / 2 + 1;
    // 一次性整体处理全部数据，保证频谱正确性
    const complex_buffer = try allocator.alloc(std.math.Complex(T), n);
    defer allocator.free(complex_buffer);
    for (0..n) |i| {
        complex_buffer[i] = std.math.Complex(T){ .re = input[i], .im = @as(T, 0.0) };
    }
    try fftInPlace(T, allocator, complex_buffer);
    convertToOutput(T, complex_buffer[0..out_len], output, magnitude);
}

const expect = std.testing.expect;
const expectApproxEqRel = std.testing.expectApproxEqRel;

fn testR2CGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;

    const size = 128;
    const input = try allocator.alloc(T, size);
    defer allocator.free(input);

    for (0..size) |i| {
        const t = @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(size));
        input[i] = @sin(2.0 * std.math.pi * @as(T, 5.0) * t) + @as(T, 0.5) * @cos(2.0 * std.math.pi * @as(T, 10.0) * t);
    }

    const out_len = size / 2 + 1;
    const output = try allocator.alloc(T, 2 * out_len);
    defer allocator.free(output);
    const magnitude = try allocator.alloc(T, out_len);
    defer allocator.free(magnitude);

    try fftR2C(T, allocator, input, output, magnitude);
    try expect(magnitude[5] > @as(T, 0.0));
    try expect(magnitude[10] > @as(T, 0.0));
}

test "Real-to-complex FFT f32" {
    try testR2CGeneric(f32);
}
test "Real-to-complex FFT f64" {
    try testR2CGeneric(f64);
}

fn testR2CSineGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tolerance: T = if (T == f32) @as(T, 1e-5) else @as(T, 1e-10);
    const N: usize = 64;
    const f: T = @as(T, 1.0);
    const input = try allocator.alloc(T, N);
    defer allocator.free(input);
    for (0..N) |i| {
        input[i] = @sin(2.0 * std.math.pi * f * @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(N)));
    }
    const out_len = N / 2 + 1;
    const output = try allocator.alloc(T, 2 * out_len);
    defer allocator.free(output);
    const magnitude = try allocator.alloc(T, out_len);
    defer allocator.free(magnitude);
    try fftR2C(T, allocator, input, output, magnitude);
    try expectApproxEqRel(magnitude[1], @as(T, @floatFromInt(N)) / @as(T, 2.0), tolerance);
}

test "fftR2C 单一正弦波 f32" {
    try testR2CSineGeneric(f32);
}
test "fftR2C 单一正弦波 f64" {
    try testR2CSineGeneric(f64);
}

fn testR2CDcGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tolerance: T = if (T == f32) @as(T, 1e-5) else @as(T, 1e-10);
    const N: usize = 32;
    const input = try allocator.alloc(T, N);
    defer allocator.free(input);
    for (0..N) |i| {
        input[i] = @as(T, 3.0);
    }
    const out_len = N / 2 + 1;
    const output = try allocator.alloc(T, 2 * out_len);
    defer allocator.free(output);
    const magnitude = try allocator.alloc(T, out_len);
    defer allocator.free(magnitude);
    try fftR2C(T, allocator, input, output, magnitude);
    try expectApproxEqRel(output[0], @as(T, @floatFromInt(N)) * @as(T, 3.0), tolerance);
    try expectApproxEqRel(output[1], @as(T, 0.0), tolerance);
}

test "fftR2C 直流分量 f32" {
    try testR2CDcGeneric(f32);
}
test "fftR2C 直流分量 f64" {
    try testR2CDcGeneric(f64);
}

fn testR2CZeroGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tolerance: T = if (T == f32) @as(T, 1e-6) else @as(T, 1e-12);
    const N: usize = 20;
    const input = try allocator.alloc(T, N);
    defer allocator.free(input);
    for (0..N) |i| {
        input[i] = @as(T, 0.0);
    }
    const out_len = N / 2 + 1;
    const output = try allocator.alloc(T, 2 * out_len);
    defer allocator.free(output);
    const magnitude = try allocator.alloc(T, out_len);
    defer allocator.free(magnitude);
    try fftR2C(T, allocator, input, output, magnitude);
    for (0..2 * out_len) |i| {
        try expectApproxEqRel(output[i], @as(T, 0.0), tolerance);
    }
}

test "fftR2C 全零输入 f32" {
    try testR2CZeroGeneric(f32);
}
test "fftR2C 全零输入 f64" {
    try testR2CZeroGeneric(f64);
}

fn testR2CEdgeGeneric(comptime T: type) !void {
    const tolerance: T = if (T == f32) @as(T, 1e-6) else @as(T, 1e-12);
    var empty: [0]T = undefined;
    var out_empty: [2]T = undefined;
    var mag_empty: [1]T = undefined;
    try fftR2C(T, std.heap.page_allocator, empty[0..], out_empty[0..], mag_empty[0..]); // 应不报错

    var one = [_]T{@as(T, 7.0)};
    var out_one: [2]T = undefined;
    var mag_one: [1]T = undefined;
    try fftR2C(T, std.heap.page_allocator, one[0..], out_one[0..], mag_one[0..]);
    try expectApproxEqRel(out_one[0], @as(T, 7.0), tolerance);
    try expectApproxEqRel(out_one[1], @as(T, 0.0), tolerance);
    try expectApproxEqRel(mag_one[0], @as(T, 7.0), tolerance);
}

test "R2C FFT edge cases f32" {
    try testR2CEdgeGeneric(f32);
}
test "R2C FFT edge cases f64" {
    try testR2CEdgeGeneric(f64);
}
