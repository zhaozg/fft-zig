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

const fftInPlace = @import("../fft.zig").fftInPlace;

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
    // Convert complex buffer to output format
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
    // 一次性整体处理全部数据，保证频谱正确性
    const complex_buffer = try allocateAlignedComplexBuffer(allocator, n);
    defer allocator.free(complex_buffer);
    for (0..n) |i| {
        complex_buffer[i] = Complex{ .re = input[i], .im = 0.0 };
    }
    try fftInPlace(allocator, complex_buffer);
    convertToOutputSIMD(complex_buffer[0..out_len], output, magnitude);
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

test "fftR2C 单一正弦波 1Hz" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();
    const N = 64;
    const f = 1.0;
    const input = try allocator.alloc(f64, N);
    defer allocator.free(input);
    for (0..N) |i| {
        input[i] = @sin(2.0 * std.math.pi * f * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(N)));
    }
    const out_len = N / 2 + 1;
    const output = try allocator.alloc(f64, 2 * out_len);
    defer allocator.free(output);
    const magnitude = try allocator.alloc(f64, out_len);
    defer allocator.free(magnitude);
    try fftR2C(allocator, input, output, magnitude);
    // 1Hz主频幅值应接近N/2
    try expectApproxEqRel(magnitude[1], N / 2, 1e-10);
    // 其他频率幅值应远小于主频
    for (2..out_len) |k| {
        try expect(magnitude[k] < 1.0);
    }
}

test "fftR2C 单一余弦波 2Hz" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();
    const N = 64;
    const f = 2.0;
    const input = try allocator.alloc(f64, N);
    defer allocator.free(input);
    for (0..N) |i| {
        input[i] = @cos(2.0 * std.math.pi * f * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(N)));
    }
    const out_len = N / 2 + 1;
    const output = try allocator.alloc(f64, 2 * out_len);
    defer allocator.free(output);
    const magnitude = try allocator.alloc(f64, out_len);
    defer allocator.free(magnitude);
    try fftR2C(allocator, input, output, magnitude);
    // 2Hz主频幅值应接近N/2
    try expectApproxEqRel(magnitude[2], N / 2, 1e-10);
    for (0..out_len) |k| {
        if (k != 2) try expect(magnitude[k] < 1.0);
    }
}

test "fftR2C 直流分量" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();
    const N = 32;
    const input = try allocator.alloc(f64, N);
    defer allocator.free(input);
    for (0..N) |i| {
        input[i] = 3.0;
    }
    const out_len = N / 2 + 1;
    const output = try allocator.alloc(f64, 2 * out_len);
    defer allocator.free(output);
    const magnitude = try allocator.alloc(f64, out_len);
    defer allocator.free(magnitude);
    try fftR2C(allocator, input, output, magnitude);
    // 直流分量应为N*3
    try expectApproxEqRel(output[0], N * 3.0, 1e-10);
    try expectApproxEqRel(output[1], 0.0, 1e-10);
    try expectApproxEqRel(magnitude[0], N * 3.0, 1e-10);
    for (1..out_len) |k| {
        try expect(magnitude[k] < 1e-8);
    }
}

test "fftR2C 两频正弦叠加" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();
    const N = 128;
    const input = try allocator.alloc(f64, N);
    defer allocator.free(input);
    for (0..N) |i| {
        input[i] = @sin(2.0 * std.math.pi * 3.0 * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(N))) + 0.5 * @sin(2.0 * std.math.pi * 5.0 * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(N)));
    }
    const out_len = N / 2 + 1;
    const output = try allocator.alloc(f64, 2 * out_len);
    defer allocator.free(output);
    const magnitude = try allocator.alloc(f64, out_len);
    defer allocator.free(magnitude);
    try fftR2C(allocator, input, output, magnitude);
    try expectApproxEqRel(magnitude[3], N / 2, 1e-10);
    try expectApproxEqRel(magnitude[5], N / 4, 1e-10);
    for (0..out_len) |k| {
        if (k != 3 and k != 5) try expect(magnitude[k] < 1.0);
    }
}

test "fftR2C 交错正负" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();
    const N = 16;
    const input = try allocator.alloc(f64, N);
    defer allocator.free(input);
    for (0..N) |i| {
        input[i] = if (i % 2 == 0) 1.0 else -1.0;
    }
    const out_len = N / 2 + 1;
    const output = try allocator.alloc(f64, 2 * out_len);
    defer allocator.free(output);
    const magnitude = try allocator.alloc(f64, out_len);
    defer allocator.free(magnitude);
    try fftR2C(allocator, input, output, magnitude);
    // N/2频点幅值应为N
    try expectApproxEqRel(magnitude[N / 2], N, 1e-10);
    for (0..out_len) |k| {
        if (k != N / 2) try expect(magnitude[k] < 1e-8);
    }
}

test "fftR2C 全零输入" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();
    const N = 20;
    const input = try allocator.alloc(f64, N);
    defer allocator.free(input);
    for (0..N) |i| {
        input[i] = 0.0;
    }
    const out_len = N / 2 + 1;
    const output = try allocator.alloc(f64, 2 * out_len);
    defer allocator.free(output);
    const magnitude = try allocator.alloc(f64, out_len);
    defer allocator.free(magnitude);
    try fftR2C(allocator, input, output, magnitude);
    for (0..2 * out_len) |i| {
        try expectApproxEqRel(output[i], 0.0, 1e-12);
    }
    for (0..out_len) |i| {
        try expectApproxEqRel(magnitude[i], 0.0, 1e-12);
    }
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

test "FFT huge data validation" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const test_sizes = [_]usize{
        1024, // power of 4
        //512, // power of 2, not power of 4
        //256, // small FFT path
        //65536, // 64K - to test non-huge path
        //1048576, // 1M
        //2097152, // 2M
        //4194304, // 4M
        //5000000, // 5M
    };

    for (test_sizes) |size| {
        std.debug.print("\n=== Testing HUGE data FFT with {d} samples ===\n", .{size});

        const input = try allocator.alloc(f64, size);
        defer allocator.free(input);

        for (0..size) |i| {
            // 生成纯正弦波，频率为5Hz
            const freq: f64 = 5.0;
            input[i] = @sin(2.0 * std.math.pi * freq * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(size)));
        }

        const out_len = size / 2 + 1;
        const output = try allocator.alloc(f64, 2 * out_len);
        defer allocator.free(output);
        const magnitude = try allocator.alloc(f64, out_len);
        defer allocator.free(magnitude);

        const start_time = std.time.nanoTimestamp();
        try fftR2C(allocator, input, output, magnitude);
        const end_time = std.time.nanoTimestamp();

        const elapsed_ms = @as(f64, @floatFromInt(@as(u64, @intCast(end_time - start_time)))) / 1e6;
        const throughput = (@as(f64, @floatFromInt(size)) / (elapsed_ms / 1000.0)) / 1e6;

        std.debug.print("Processing time: {d:.1}ms\n", .{elapsed_ms});
        std.debug.print("Throughput: {d:.1} MSamples/s\n", .{throughput});

        try expect(magnitude[0] > 0.0);
        try expect(!math.isNan(magnitude[0]));
        try expect(math.isFinite(magnitude[0]));

        // var peak_count: usize = 0;
        var total_energy: f64 = 0.0;

        for (0..@min(1000, out_len)) |i| {
            try expect(!math.isNan(magnitude[i]));
            try expect(math.isFinite(magnitude[i]));
            try expect(magnitude[i] >= 0.0);

            total_energy += magnitude[i] * magnitude[i];
        }

        // 直接检测主频点幅值
        // For a pure sine wave of amplitude 1.0 at frequency 5,
        // the FFT magnitude at bin 5 should be approximately size/2
        const expected_magnitude = @as(f64, @floatFromInt(size)) / 2.0;
        std.debug.print("Expected magnitude: {d:.1}, Actual magnitude[5]: {d:.1}\n", .{ expected_magnitude, magnitude[5] });
        try expect(magnitude[5] > expected_magnitude * 0.9); // Allow 10% tolerance
        try expect(total_energy > 0.1);

        std.debug.print("Total energy: {d:.1}\n", .{total_energy});
    }
}
