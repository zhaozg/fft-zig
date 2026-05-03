//! Performance benchmarking tests for FFT algorithms
//! 性能测试与基准测试

const std = @import("std");
const fftInPlace = @import("base.zig").fftInPlaceBase;
const fft_radix2 = @import("fft_radix2.zig");
const fft_radix4 = @import("fft_radix4.zig");
const fft_utils = @import("utils.zig");

const expect = std.testing.expect;

// 性能基准测试：比较不同大小的FFT性能
fn testBenchmarkGeneric(comptime T: type) !void {
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    const clock = std.Io.Clock.awake;

    const test_sizes = [_]usize{ 64, 128, 256, 512, 1024, 2048, 4096, 8192 };

    std.debug.print("\n=== FFT Performance Benchmark ({s}) ===\n", .{@typeName(T)});
    for (test_sizes) |size| {
        var data = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(data);

        var prng = std.Random.DefaultPrng.init(42);
        const random = prng.random();
        for (0..size) |i| {
            data[i] = std.math.Complex(T){
                .re = random.float(T) * @as(T, 2.0) - @as(T, 1.0),
                .im = @as(T, 0.0),
            };
        }

        const start_time = std.Io.Clock.now(clock, io).toNanoseconds();
        try fftInPlace(T, allocator, data);
        const end_time = std.Io.Clock.now(clock, io).toNanoseconds();

        const elapsed_s = @as(f128, @floatFromInt(@as(i96, @intCast(end_time - start_time)))) / std.time.us_per_s;
        const throughput = (@as(f128, @floatFromInt(size)) / (elapsed_s) / 1e6);

        std.debug.print("Size {d:6}: {d:9.3} μs, {d:9.3} MSamples/s\n", .{
            size,
            elapsed_s,
            throughput,
        });
    }
    std.debug.print("\n", .{});
}

test "FFT performance benchmark: Various sizes f32" {
    try testBenchmarkGeneric(f32);
}
test "FFT performance benchmark: Various sizes f64" {
    try testBenchmarkGeneric(f64);
}

// 比较Radix-2和Radix-4性能
fn testRadixCompareGeneric(comptime T: type) !void {
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    const clock = std.Io.Clock.awake;
    // TODO: threshold 缩小阀值
    const threshold = if (T == f64) 1e-5 else 1e-2;

    const test_sizes = [_]usize{ 64, 256, 1024, 4096 };

    std.debug.print("\n=== Radix-2 vs Radix-4 Performance ({s}) ===\n", .{@typeName(T)});
    for (test_sizes) |size| {
        try expect(fft_utils.isPowerOfFour(size));

        var data2 = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(data2);
        var data4 = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(data4);

        var prng = std.Random.DefaultPrng.init(12345);
        const random = prng.random();
        for (0..size) |i| {
            const val = std.math.Complex(T){
                .re = random.float(T) * @as(T, 2.0) - @as(T, 1.0),
                .im = @as(T, 0.0),
            };
            data2[i] = val;
            data4[i] = val;
        }

        var start_time = std.Io.Clock.now(clock, io).toNanoseconds();
        try fft_radix2.fftRadix2SIMD(T, data2);
        var end_time = std.Io.Clock.now(clock, io).toNanoseconds();
        const time2_ns = end_time - start_time;

        start_time = std.Io.Clock.now(clock, io).toNanoseconds();
        try fft_radix4.fftRadix4SIMD(T, data4);
        end_time = std.Io.Clock.now(clock, io).toNanoseconds();
        const time4_ns = end_time - start_time;

        const speedup = @as(f128, @floatFromInt(time2_ns)) / @as(f128, @floatFromInt(time4_ns));

        std.debug.print("Size {d:6}: Radix-2={d:9}ns, Radix-4={d:9}ns, Speedup={d:9.3}x\n", .{
            size,
            time2_ns,
            time4_ns,
            speedup,
        });

        for (0..size) |i| {
            const diff_re = @abs(data2[i].re - data4[i].re);
            const diff_im = @abs(data2[i].im - data4[i].im);
            try expect(diff_re < @as(T, threshold));
            try expect(diff_im < @as(T, threshold));
        }
    }
    std.debug.print("\n", .{});
}

test "FFT performance: Radix-2 vs Radix-4 f32" {
    try testRadixCompareGeneric(f32);
}
test "FFT performance: Radix-2 vs Radix-4 f64" {
    try testRadixCompareGeneric(f64);
}

// 测试大数据处理性能
fn testLargeDataGeneric(comptime T: type) !void {
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    const clock = std.Io.Clock.awake;

    const large_sizes = [_]usize{ 16384, 65536, 262144 };

    std.debug.print("\n=== Large Data FFT Performance ({s}) ===\n", .{@typeName(T)});
    for (large_sizes) |size| {
        var data = try allocator.alloc(std.math.Complex(T), size);
        defer allocator.free(data);

        var prng = std.Random.DefaultPrng.init(99999);
        const random = prng.random();
        for (0..size) |i| {
            data[i] = std.math.Complex(T){
                .re = random.float(T) * @as(T, 2.0) - @as(T, 1.0),
                .im = @as(T, 0.0),
            };
        }

        const start_time = std.Io.Clock.now(clock, io).toNanoseconds();
        try fftInPlace(T, allocator, data);
        const end_time = std.Io.Clock.now(clock, io).toNanoseconds();

        const elapsed_s = @as(f128, @floatFromInt(@as(i96, @intCast(end_time - start_time)))) / std.time.us_per_s;
        const throughput = (@as(f128, @floatFromInt(size)) / (elapsed_s) / 1e6);

        std.debug.print("Size {d:6}: {d:9.3} s, {d:9.3} MSamples/s\n", .{
            size,
            elapsed_s,
            throughput,
        });

        var total_energy: T = @as(T, 0.0);
        for (data) |v| {
            total_energy += v.re * v.re + v.im * v.im;
        }
        try expect(total_energy > @as(T, 0.0));
    }
    std.debug.print("\n", .{});
}

test "FFT performance: Large data handling f32" {
    try testLargeDataGeneric(f32);
}
test "FFT performance: Large data handling f64" {
    try testLargeDataGeneric(f64);
}

// 测试内存分配效率
fn testMemoryGeneric(comptime T: type) !void {
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    const clock = std.Io.Clock.awake;

    const size: usize = 4096;
    const iterations: usize = 100;

    std.debug.print("\n=== Memory Allocation Performance ({s}) ===\n", .{@typeName(T)});

    var total_time_ns: u64 = 0;
    for (0..iterations) |_| {
        const start_time = std.Io.Clock.now(clock, io).toNanoseconds();

        var data = try allocator.alloc(std.math.Complex(T), size);
        for (0..size) |i| {
            data[i] = std.math.Complex(T){ .re = @as(T, @floatFromInt(i)), .im = @as(T, 0.0) };
        }
        try fftInPlace(T, allocator, data);
        allocator.free(data);

        const end_time = std.Io.Clock.now(clock, io).toNanoseconds();
        total_time_ns += @as(u64, @intCast(end_time - start_time));
    }

    const avg_time_us = @as(f128, @floatFromInt(total_time_ns)) / @as(f128, @floatFromInt(iterations));
    std.debug.print("Average time per iteration (size {d:6}): {d:9.2} μs\n", .{ size, avg_time_us });
    std.debug.print("Total iterations: {d:6}\n", .{iterations});
    std.debug.print("\n", .{});
}

test "FFT performance: Memory allocation patterns f32" {
    try testMemoryGeneric(f32);
}
test "FFT performance: Memory allocation patterns f64" {
    try testMemoryGeneric(f64);
}
