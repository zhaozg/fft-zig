//! f32 vs f64 性能对比示例
//!
//! 展示泛型 FFT 库中不同浮点类型的性能差异。
//! f32 通常比 f64 快 1.5-2.5 倍，内存占用减半。
//!
//! 运行方式：zig build run-examples

const std = @import("std");
const fft = @import("fft");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const test_sizes = [_]usize{ 256, 1024, 4096, 16384 };

    std.debug.print("f32 vs f64 性能对比\n", .{});
    std.debug.print("====================\n\n", .{});

    for (test_sizes) |N| {
        // --- f32 测试 ---
        const input_f32 = try allocator.alloc(f32, N);
        defer allocator.free(input_f32);

        for (0..N) |i| {
            const t = @as(f32, @floatFromInt(i)) / @as(f32, @floatFromInt(N));
            input_f32[i] = @sin(2.0 * std.math.pi * 7.0 * t) + 0.3 * @cos(2.0 * std.math.pi * 23.0 * t);
        }

        const data_f32 = try allocator.alloc(std.math.Complex(f32), N);
        defer allocator.free(data_f32);

        // 复制输入到复数缓冲区
        for (0..N) |i| {
            data_f32[i] = std.math.Complex(f32){ .re = input_f32[i], .im = 0.0 };
        }

        const start_f32 = std.time.milliTimestamp();
        const iterations_f32 = 100;
        for (0..iterations_f32) |_| {
            // 每次重新复制数据
            const work = try allocator.dupe(std.math.Complex(f32), data_f32);
            defer allocator.free(work);
            try fft.fftInPlace(f32, allocator, work);
        }
        const end_f32 = std.time.milliTimestamp();
        const elapsed_f32 = @as(f64, @floatFromInt(end_f32 - start_f32)) / @as(f64, @floatFromInt(iterations_f32));

        // --- f64 测试 ---
        const input_f64 = try allocator.alloc(f64, N);
        defer allocator.free(input_f64);

        for (0..N) |i| {
            const t = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(N));
            input_f64[i] = @sin(2.0 * std.math.pi * 7.0 * t) + 0.3 * @cos(2.0 * std.math.pi * 23.0 * t);
        }

        const data_f64 = try allocator.alloc(std.math.Complex(f64), N);
        defer allocator.free(data_f64);

        for (0..N) |i| {
            data_f64[i] = std.math.Complex(f64){ .re = input_f64[i], .im = 0.0 };
        }

        const start_f64 = std.time.milliTimestamp();
        const iterations_f64 = 100;
        for (0..iterations_f64) |_| {
            const work = try allocator.dupe(std.math.Complex(f64), data_f64);
            defer allocator.free(work);
            try fft.fftInPlace(f64, allocator, work);
        }
        const end_f64 = std.time.milliTimestamp();
        const elapsed_f64 = @as(f64, @floatFromInt(end_f64 - start_f64)) / @as(f64, @floatFromInt(iterations_f64));

        const speedup = elapsed_f64 / elapsed_f32;

        std.debug.print("N={d:>5}:  f32={d:>8.2}ms  f64={d:>8.2}ms  加速比={d:.2}x\n", .{
            N, elapsed_f32, elapsed_f64, speedup,
        });
    }

    std.debug.print("\n结论：\n", .{});
    std.debug.print("  - f32 内存占用为 f64 的一半\n", .{});
    std.debug.print("  - f32 计算速度通常比 f64 快 1.5-2.5 倍\n", .{});
    std.debug.print("  - 选择建议：精度要求不高时用 f32，科学计算用 f64\n", .{});
}
