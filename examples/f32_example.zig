//! f32 FFT 示例 - 展示如何使用泛型 FFT 库处理单精度浮点数据
//!
//! 使用 f32 的优势：
//! - 内存占用减半（相比 f64）
//! - 计算速度更快（SIMD 可处理更多元素）
//! - 适合音频、图像等对精度要求不高的场景
//!
//! 运行方式：zig build run-examples

const std = @import("std");
const fft = @import("fft");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const N = 1024; // FFT 点数

    // 准备输入数据：一个 7Hz 的正弦波 + 0.3 倍 23Hz 的余弦波
    const input = try allocator.alloc(f32, N);
    defer allocator.free(input);

    for (0..N) |i| {
        const t = @as(f32, @floatFromInt(i)) / @as(f32, @floatFromInt(N));
        input[i] = @sin(2.0 * std.math.pi * 7.0 * t) + 0.3 * @cos(2.0 * std.math.pi * 23.0 * t);
    }

    // 分配输出缓冲区
    const output = try allocator.alloc(std.math.Complex(f32), N);
    defer allocator.free(output);

    // 调用泛型 fft 函数，指定类型为 f32
    try fft.fft(f32, allocator, input, output);

    // 查找峰值频率
    var max_magnitude: f32 = 0.0;
    var peak_freq: usize = 0;
    for (1..N / 2) |i| {
        const magnitude = @sqrt(output[i].re * output[i].re + output[i].im * output[i].im);
        if (magnitude > max_magnitude) {
            max_magnitude = magnitude;
            peak_freq = i;
        }
    }

    std.debug.print("f32 FFT 结果：\n", .{});
    std.debug.print("  输入信号：sin(2π·7·t) + 0.3·cos(2π·23·t)\n", .{});
    std.debug.print("  采样点数：{d}\n", .{N});
    std.debug.print("  峰值频率索引：{d}（预期：7 和 23）\n", .{peak_freq});
    std.debug.print("  峰值幅度：{d:.4}\n", .{max_magnitude});

    // 验证正确性
    const expected_freq = 7;
    if (peak_freq >= expected_freq - 2 and peak_freq <= expected_freq + 2) {
        std.debug.print("  ✅ 频率检测正确！\n", .{});
    } else {
        std.debug.print("  ❌ 频率检测异常！\n", .{});
    }
}
