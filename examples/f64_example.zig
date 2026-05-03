//! f64 FFT 示例 - 展示如何使用泛型 FFT 库处理双精度浮点数据
//!
//! 使用 f64 的优势：
//! - 更高精度（约 15-17 位有效数字）
//! - 适合科学计算、信号处理等对精度要求高的场景
//! - 与旧版 API 完全兼容（fft_f64 别名）
//!
//! 运行方式：zig build run-examples

const std = @import("std");
const fft = @import("fft");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const N = 2048; // FFT 点数

    // 准备输入数据：一个复合信号
    const input = try allocator.alloc(f64, N);
    defer allocator.free(input);

    for (0..N) |i| {
        const t = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(N));
        // 50Hz 主频 + 120Hz 谐波 + 噪声
        input[i] = 0.7 * @sin(2.0 * std.math.pi * 50.0 * t) +
                   0.3 * @sin(2.0 * std.math.pi * 120.0 * t) +
                   0.1 * @sin(2.0 * std.math.pi * 200.0 * t);
    }

    // 分配输出缓冲区
    const output = try allocator.alloc(std.math.Complex(f64), N);
    defer allocator.free(output);

    // 方式一：使用泛型 API（推荐）
    try fft.fft(f64, allocator, input, output);

    // 方式二：使用向后兼容别名（效果相同）
    // try fft.fft_f64(allocator, input, output);

    // 分析频谱
    std.debug.print("f64 FFT 频谱分析：\n", .{});
    std.debug.print("  输入信号：0.7·sin(2π·50·t) + 0.3·sin(2π·120·t) + 0.1·sin(2π·200·t)\n", .{});
    std.debug.print("  采样点数：{d}\n", .{N});
    std.debug.print("  频率分辨率：{d:.2} Hz\n", .{@as(f64, 1.0)});

    // 查找所有显著频率分量
    std.debug.print("\n  检测到的频率分量：\n", .{});
    for (1..N / 2) |i| {
        const magnitude = @sqrt(output[i].re * output[i].re + output[i].im * output[i].im);
        if (magnitude > 10.0) { // 阈值筛选
            std.debug.print("    频率索引 {d:>4}：幅度 {d:.2}\n", .{ i, magnitude });
        }
    }

    // 验证 Parseval 定理（能量守恒）
    var energy_time: f64 = 0.0;
    for (input) |x| {
        energy_time += x * x;
    }
    var energy_freq: f64 = 0.0;
    for (output) |X| {
        energy_freq += X.re * X.re + X.im * X.im;
    }
    energy_freq /= @as(f64, @floatFromInt(N));

    std.debug.print("\n  能量守恒验证（Parseval 定理）：\n", .{});
    std.debug.print("    时域能量：{d:.10}\n", .{energy_time});
    std.debug.print("    频域能量：{d:.10}\n", .{energy_freq});
    std.debug.print("    相对误差：{d:.2e}\n", .{@abs(energy_time - energy_freq) / energy_time});
}
