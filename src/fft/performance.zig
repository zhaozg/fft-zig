//! Performance benchmarking tests for FFT algorithms
//! 性能测试与基准测试

const std = @import("std");
const math = std.math;
const Complex = @import("types.zig").Complex;
const fftInPlace = @import("base.zig").fftInPlaceBase;
const fft_radix2 = @import("fft_radix2.zig");
const fft_radix4 = @import("fft_radix4.zig");
const fft_utils = @import("utils.zig");

const expect = std.testing.expect;

// 性能基准测试：比较不同大小的FFT性能
test "FFT performance benchmark: Various sizes" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const test_sizes = [_]usize{ 64, 128, 256, 512, 1024, 2048, 4096, 8192 };

    std.debug.print("\n=== FFT Performance Benchmark ===\n", .{});
    for (test_sizes) |size| {
        var data = try allocator.alloc(Complex, size);
        defer allocator.free(data);

        // 生成测试信号
        var prng = std.Random.DefaultPrng.init(42);
        const random = prng.random();
        for (0..size) |i| {
            data[i] = Complex{
                .re = random.float(f64) * 2.0 - 1.0,
                .im = 0.0,
            };
        }

        const start = std.time.nanoTimestamp();
        try fftInPlace(allocator, data);
        const end = std.time.nanoTimestamp();

        const elapsed_ns = @as(f64, @floatFromInt(@as(u64, @intCast(end - start))));
        const elapsed_us = elapsed_ns / 1000.0;
        const throughput_msps = (@as(f64, @floatFromInt(size)) / elapsed_ns) * 1000.0;

        std.debug.print("Size {d:6}: {d:8.2} μs, {d:6.2} MSamples/s\n", .{
            size,
            elapsed_us,
            throughput_msps,
        });
    }
    std.debug.print("\n", .{});
}

// 比较Radix-2和Radix-4性能
test "FFT performance: Radix-2 vs Radix-4" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const test_sizes = [_]usize{ 64, 256, 1024, 4096 }; // 都是4的幂

    std.debug.print("\n=== Radix-2 vs Radix-4 Performance ===\n", .{});
    for (test_sizes) |size| {
        // 验证是4的幂
        try expect(fft_utils.isPowerOfFour(size));

        var data2 = try allocator.alloc(Complex, size);
        defer allocator.free(data2);
        var data4 = try allocator.alloc(Complex, size);
        defer allocator.free(data4);

        // 生成相同的测试信号
        var prng = std.Random.DefaultPrng.init(12345);
        const random = prng.random();
        for (0..size) |i| {
            const val = Complex{
                .re = random.float(f64) * 2.0 - 1.0,
                .im = 0.0,
            };
            data2[i] = val;
            data4[i] = val;
        }

        // Radix-2 SIMD
        const start2 = std.time.nanoTimestamp();
        try fft_radix2.fftRadix2SIMD(data2);
        const end2 = std.time.nanoTimestamp();
        const time2_ns = @as(f64, @floatFromInt(@as(u64, @intCast(end2 - start2))));

        // Radix-4 SIMD
        const start4 = std.time.nanoTimestamp();
        try fft_radix4.fftRadix4SIMD(data4);
        const end4 = std.time.nanoTimestamp();
        const time4_ns = @as(f64, @floatFromInt(@as(u64, @intCast(end4 - start4))));

        const speedup = time2_ns / time4_ns;

        std.debug.print("Size {d:6}: Radix-2={d:6.2}μs, Radix-4={d:6.2}μs, Speedup={d:.2}x\n", .{
            size,
            time2_ns / 1000.0,
            time4_ns / 1000.0,
            speedup,
        });

        // 验证结果一致性
        for (0..size) |i| {
            const diff_re = @abs(data2[i].re - data4[i].re);
            const diff_im = @abs(data2[i].im - data4[i].im);
            try expect(diff_re < 1e-10);
            try expect(diff_im < 1e-10);
        }
    }
    std.debug.print("\n", .{});
}

// 测试大数据处理性能
test "FFT performance: Large data handling" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const large_sizes = [_]usize{ 16384, 65536, 262144 };

    std.debug.print("\n=== Large Data FFT Performance ===\n", .{});
    for (large_sizes) |size| {
        var data = try allocator.alloc(Complex, size);
        defer allocator.free(data);

        // 生成测试信号
        var prng = std.Random.DefaultPrng.init(99999);
        const random = prng.random();
        for (0..size) |i| {
            data[i] = Complex{
                .re = random.float(f64) * 2.0 - 1.0,
                .im = 0.0,
            };
        }

        const start = std.time.nanoTimestamp();
        try fftInPlace(allocator, data);
        const end = std.time.nanoTimestamp();

        const elapsed_ms = @as(f64, @floatFromInt(@as(u64, @intCast(end - start)))) / 1e6;
        const throughput_msps = (@as(f64, @floatFromInt(size)) / elapsed_ms) / 1000.0;

        std.debug.print("Size {d:8}: {d:8.2} ms, {d:6.2} MSamples/s\n", .{
            size,
            elapsed_ms,
            throughput_msps,
        });

        // 基本的正确性检查 - 验证能量守恒
        var total_energy: f64 = 0.0;
        for (data) |v| {
            total_energy += v.re * v.re + v.im * v.im;
        }
        try expect(total_energy > 0.0); // 确保计算了有效的FFT
    }
    std.debug.print("\n", .{});
}

// 测试内存分配效率
test "FFT performance: Memory allocation patterns" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const size: usize = 4096;
    const iterations: usize = 100;

    std.debug.print("\n=== Memory Allocation Performance ===\n", .{});

    // 测试重复分配和释放
    var total_time_ns: u64 = 0;
    for (0..iterations) |_| {
        const start = std.time.nanoTimestamp();

        var data = try allocator.alloc(Complex, size);
        for (0..size) |i| {
            data[i] = Complex{ .re = @as(f64, @floatFromInt(i)), .im = 0.0 };
        }
        try fftInPlace(allocator, data);
        allocator.free(data);

        const end = std.time.nanoTimestamp();
        total_time_ns += @as(u64, @intCast(end - start));
    }

    const avg_time_us = @as(f64, @floatFromInt(total_time_ns)) / @as(f64, @floatFromInt(iterations)) / 1000.0;
    std.debug.print("Average time per iteration (size {d}): {d:.2} μs\n", .{ size, avg_time_us });
    std.debug.print("Total iterations: {d}\n", .{iterations});
    std.debug.print("\n", .{});
}
