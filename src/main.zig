// filepath: src/main.zig
const std = @import("std");
const fft = @import("fft.zig");

fn spectralTest(allocator: std.mem.Allocator, bits: []const f64) !f64 {
    // 离散傅里叶检验（频谱检验）- GMT0005
    const n = bits.len;

    // 分配FFT输出缓冲区
    const out_len = n / 2 + 1;
    const fft_out = try allocator.alloc(f64, 2 * out_len);
    defer allocator.free(fft_out);
    const fft_m = try allocator.alloc(f64, out_len);
    defer allocator.free(fft_m);

    // 执行FFT
    try fft.fftR2C(allocator, bits, fft_out, fft_m);

    // 计算峰值阈值
    const threshold = @sqrt(@log(1.0 / 0.05) * @as(f64, @floatFromInt(n)));

    // 统计超过阈值的峰值数量
    var peaks: usize = 0;
    for (fft_m[0 .. n / 2]) |magnitude| {
        if (magnitude < threshold) {
            peaks += 1;
        }
    }

    const n0 = @as(f64, @floatFromInt(peaks));
    const n_half = @as(f64, @floatFromInt(n)) / 2.0;
    const expected = 0.95 * n_half;

    const d = (n0 - expected) / @sqrt(n_half * 0.95 * 0.05);

    return erfc(@abs(d) / @sqrt(2.0));
}

// 数学辅助函数
fn erfc(x: f64) f64 {
    // 互补误差函数的近似实现
    const a1: f64 = 0.254829592;
    const a2: f64 = -0.284496736;
    const a3: f64 = 1.421413741;
    const a4: f64 = -1.453152027;
    const a5: f64 = 1.061405429;
    const p: f64 = 0.3275911;

    const sign: f64 = if (x < 0) -1.0 else 1.0;
    const abs_x = @abs(x);

    const t = 1.0 / (1.0 + p * abs_x);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * @exp(-abs_x * abs_x);

    return if (sign > 0) y else 2.0 - y;
}

pub fn main() !void {
    const allocator = std.heap.page_allocator;
    var args = std.process.args();
    defer args.deinit();

    _ = args.next(); // 跳过程序名
    const filename = args.next() orelse {
        std.debug.print("用法: fft_bin <input_file>\n", .{});
        std.debug.print("功能: 对二进制文件进行FFT频谱分析和离散傅里叶检验\n", .{});
        return error.InvalidArgs;
    };

    // 打开并读取文件
    var file = try std.fs.cwd().openFile(filename, .{});
    defer file.close();

    const file_size = try file.getEndPos();
    std.debug.print("文件大小: {d} 字节, {d} 比特\n", .{ file_size, file_size * 8 });

    const nbits = file_size * 8;
    const buf = try allocator.alloc(f64, nbits);
    defer allocator.free(buf);

    var byte: u8 = undefined;
    var idx: usize = 0;

    while (file.read(std.mem.asBytes(&byte)) catch 0 != 0) {
        for (0..8) |i| {
            const bit = (byte >> @as(u3, @intCast(i))) & 1;
            buf[idx] = if (bit == 1) 1.0 else -1.0;
            idx += 1;
        }
    }

    std.debug.print("读取完成，共 {d} 比特\n", .{idx});

    // 执行GMT0005离散傅里叶检验
    std.debug.print("\n=== GMT0005 离散傅里叶检验 ===\n", .{});
    const spectral_p_value = try spectralTest(allocator, buf);
    std.debug.print("离散傅里叶检验 P值: {d:.6}\n", .{spectral_p_value});

    const threshold: f64 = 0.01; // 显著性水平
    if (spectral_p_value >= threshold) {
        std.debug.print("检验结果: ✓ 通过 (P值 >= {d:.2})\n", .{threshold});
    } else {
        std.debug.print("检验结果: ✗ 未通过 (P值 < {d:.2})\n", .{threshold});
    }

    // 执行FFT频谱分析
    std.debug.print("\n=== FFT频谱分析 ===\n", .{});

    const out_len = nbits / 2 + 1;
    const fft_out = try allocator.alloc(f64, 2 * out_len);
    defer allocator.free(fft_out);
    const fft_m = try allocator.alloc(f64, out_len);
    defer allocator.free(fft_m);

    const start_time = std.time.nanoTimestamp();
    try fft.fftR2C(allocator, buf, fft_out, fft_m);
    const end_time = std.time.nanoTimestamp();

    const elapsed_ms = @as(f64, @floatFromInt(@as(u64, @intCast(end_time - start_time)))) / 1e6;
    std.debug.print("FFT处理时间: {d:.2}ms\n", .{elapsed_ms});

    // 打印前10个幅值和统计信息
    std.debug.print("\n前10个频率分量幅值:\n", .{});
    for (0..@min(10, out_len)) |i| {
        std.debug.print("  频率 {d}: {d:.6}\n", .{ i, fft_m[i] });
    }

    // 计算频谱统计
    var max_magnitude: f64 = 0.0;
    var max_freq: usize = 0;
    var total_power: f64 = 0.0;

    for (0..out_len) |i| {
        total_power += fft_m[i] * fft_m[i];
        if (fft_m[i] > max_magnitude) {
            max_magnitude = fft_m[i];
            max_freq = i;
        }
    }

    std.debug.print("\n频谱统计:\n", .{});
    std.debug.print("  总功率: {d:.2}\n", .{total_power});
    std.debug.print("  最大幅值: {d:.6} (频率索引 {d})\n", .{ max_magnitude, max_freq });
    std.debug.print("  平均功率: {d:.6}\n", .{total_power / @as(f64, @floatFromInt(out_len))});
}
