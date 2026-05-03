const std = @import("std");
const fft = @import("fft.zig");

const expectApproxEqRel = std.testing.expectApproxEqRel;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;
const expect = std.testing.expect;

fn testFFTPerformanceGeneric(comptime T: type) !void {
    const io = std.testing.io;
    const allocator = std.testing.allocator;

    const test_sizes = [_]usize{ 256, 1024, 4096 };
    for (test_sizes) |size| {
        std.debug.print("Benchmarking size {d}...\n", .{size});

        var input = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(input);

        for (0..size) |i| {
            const t = @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(size));
            input[i] = std.math.Complex(T){ .re = @sin(2.0 * std.math.pi * @as(T, 7.0) * t) + @as(T, 0.3) * @cos(2.0 * std.math.pi * @as(T, 23.0) * t), .im = @as(T, 0.0) };
        }

        const data = try allocator.dupe(std.math.Complex(T), input);
        defer allocator.free(data);

        const clock = std.Io.Clock.awake;
        const start_time = std.Io.Clock.now(clock, io).toNanoseconds();
        try fft.fftInPlace(T, allocator, data);
        const end_time = std.Io.Clock.now(clock, io).toNanoseconds();

        const elapsed_s = @as(f128, @floatFromInt(@as(i96, @intCast(end_time - start_time)))) / std.time.us_per_s;
        const throughput = (@as(f128, @floatFromInt(size)) / (elapsed_s) / 1e6);

        std.debug.print("  Size {d}: {d:.2}ms, {d:.1} MSamples/s\n", .{ size, elapsed_s, throughput });

        // Validate correctness - find dominant frequency
        var max_magnitude: T = @as(T, 0.0);
        var peak_freq: usize = 0;
        for (1..size / 2) |i| {
            const magnitude = @sqrt(data[i].re * data[i].re + data[i].im * data[i].im);
            if (magnitude > max_magnitude) {
                max_magnitude = magnitude;
                peak_freq = i;
            }
        }
        const expected_freq = (7 * size) / size;
        try expect(peak_freq >= expected_freq - 2 and peak_freq <= expected_freq + 2);
    }
}

test "FFT performance and large data f32" {
    try testFFTPerformanceGeneric(f32);
}
test "FFT performance and large data f64" {
    try testFFTPerformanceGeneric(f64);
}

test {
    _ = @import("fft/base.zig");
    _ = @import("fft/fft_r2c.zig");
    _ = @import("fft/ifft.zig");
    _ = @import("fft/utils.zig");
    _ = @import("fft/fft_mixed.zig");
    _ = @import("fft/fft_radix2.zig");
    _ = @import("fft/twiddle.zig");
    _ = @import("fft/fft_parallel.zig");
    _ = @import("fft/fft_radix4.zig");
    _ = @import("fft/types.zig");
    _ = @import("fft/validation.zig");
    _ = @import("fft/performance.zig");
    _ = @import("fft/edge_cases.zig");
}
