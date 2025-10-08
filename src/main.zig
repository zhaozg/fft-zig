// filepath: src/main.zig
const std = @import("std");
const fft = @import("fft.zig");

// GMT0005 随机性检测指标
const RandomnessMetrics = struct {
    monobit_test: f64, // 单比特频数检验
    frequency_block_test: f64, // 块内频数检验
    runs_test: f64, // 游程检验
    longest_run_test: f64, // 块内最长游程检验
    spectral_test: f64, // 离散傅里叶检验 (DFT/FFT)

    pub fn evaluate(self: *const RandomnessMetrics) void {
        std.debug.print("\n=== GMT0005 随机性检测评价 ===\n", .{});
        std.debug.print("1. 单比特频数检验 P值: {d:.6}\n", .{self.monobit_test});
        std.debug.print("2. 块内频数检验 P值: {d:.6}\n", .{self.frequency_block_test});
        std.debug.print("3. 游程检验 P值: {d:.6}\n", .{self.runs_test});
        std.debug.print("4. 块内最长游程检验 P值: {d:.6}\n", .{self.longest_run_test});
        std.debug.print("5. 离散傅里叶检验 P值: {d:.6}\n", .{self.spectral_test});

        const threshold: f64 = 0.01; // 显著性水平
        var passed: usize = 0;
        const total: usize = 5;

        if (self.monobit_test >= threshold) passed += 1;
        if (self.frequency_block_test >= threshold) passed += 1;
        if (self.runs_test >= threshold) passed += 1;
        if (self.longest_run_test >= threshold) passed += 1;
        if (self.spectral_test >= threshold) passed += 1;

        std.debug.print("\n通过检验: {d}/{d}\n", .{ passed, total });

        if (passed == total) {
            std.debug.print("整体评价: ✓ 通过所有随机性检验\n", .{});
        } else if (passed >= 4) {
            std.debug.print("整体评价: △ 基本通过随机性检验\n", .{});
        } else {
            std.debug.print("整体评价: ✗ 未通过随机性检验\n", .{});
        }
    }
};

fn monobitTest(bits: []const f64) f64 {
    // 单比特频数检验：统计0和1的数量是否均衡
    var sum: f64 = 0.0;
    for (bits) |bit| {
        sum += bit; // +1 for 1, -1 for 0
    }

    const n = @as(f64, @floatFromInt(bits.len));
    const s_obs = @abs(sum) / @sqrt(n);

    // 计算P值（使用互补误差函数近似）
    return erfc(s_obs / @sqrt(2.0));
}

fn frequencyBlockTest(bits: []const f64, block_size: usize) f64 {
    // 块内频数检验
    const n = bits.len;
    const num_blocks = n / block_size;

    if (num_blocks == 0) return 1.0;

    var chi_squared: f64 = 0.0;
    var i: usize = 0;
    while (i < num_blocks) : (i += 1) {
        var block_sum: f64 = 0.0;
        for (0..block_size) |j| {
            block_sum += bits[i * block_size + j];
        }
        const pi = (block_sum / @as(f64, @floatFromInt(block_size)) + 1.0) / 2.0;
        chi_squared += (pi - 0.5) * (pi - 0.5);
    }

    chi_squared *= 4.0 * @as(f64, @floatFromInt(block_size));

    // 返回近似P值
    const df = @as(f64, @floatFromInt(num_blocks));
    return igamc(df / 2.0, chi_squared / 2.0);
}

fn runsTest(bits: []const f64) f64 {
    // 游程检验：检测序列中游程的数量
    if (bits.len < 2) return 1.0;

    var runs: usize = 1;
    for (1..bits.len) |i| {
        if (bits[i] != bits[i - 1]) {
            runs += 1;
        }
    }

    // 计算比例
    var ones: usize = 0;
    for (bits) |bit| {
        if (bit > 0) ones += 1;
    }

    const n = @as(f64, @floatFromInt(bits.len));
    const pi = @as(f64, @floatFromInt(ones)) / n;

    // 前提检验：pi应该接近0.5
    if (@abs(pi - 0.5) >= 2.0 / @sqrt(n)) {
        return 0.0;
    }

    const v_obs = @as(f64, @floatFromInt(runs));
    const expected = 2.0 * n * pi * (1.0 - pi);
    const variance = 2.0 * n * pi * (1.0 - pi) * (2.0 * n * pi * (1.0 - pi) - 1.0) / (n - 1.0);

    const z = (v_obs - expected) / @sqrt(variance);

    return erfc(@abs(z) / @sqrt(2.0));
}

fn longestRunTest(bits: []const f64, block_size: usize) f64 {
    // 块内最长游程检验
    const n = bits.len;
    const num_blocks = n / block_size;

    if (num_blocks == 0) return 1.0;

    var i: usize = 0;
    var max_run_sum: f64 = 0.0;

    while (i < num_blocks) : (i += 1) {
        var max_run: usize = 0;
        var current_run: usize = 0;
        var prev_bit: f64 = 0;

        for (0..block_size) |j| {
            const bit = bits[i * block_size + j];
            if (bit == prev_bit and bit > 0) {
                current_run += 1;
            } else {
                if (current_run > max_run) max_run = current_run;
                current_run = if (bit > 0) 1 else 0;
            }
            prev_bit = bit;
        }
        if (current_run > max_run) max_run = current_run;

        max_run_sum += @as(f64, @floatFromInt(max_run));
    }

    const avg_max_run = max_run_sum / @as(f64, @floatFromInt(num_blocks));
    const expected = @log2(@as(f64, @floatFromInt(block_size)));

    // 简化的P值计算
    const deviation = @abs(avg_max_run - expected) / expected;
    return @max(0.0, 1.0 - deviation);
}

fn spectralTest(allocator: std.mem.Allocator, bits: []const f64) !f64 {
    // 离散傅里叶检验（频谱检验）
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

fn igamc(a: f64, x: f64) f64 {
    // 不完全伽马函数的互补函数近似
    // 简化实现，返回近似P值
    if (x < 0 or a <= 0) return 1.0;
    if (x == 0) return 1.0;

    // 使用近似公式
    const lambda = x / a;
    if (lambda < 1) {
        return 1.0 - @exp(-x) * std.math.pow(f64, x, a) / gamma(a + 1);
    } else {
        return @exp(-x) * std.math.pow(f64, x, a) / gamma(a);
    }
}

fn gamma(x: f64) f64 {
    // Stirling近似
    if (x < 0.5) return std.math.pi / (@sin(std.math.pi * x) * gamma(1.0 - x));
    const z = x - 1.0;
    const sqrt_2pi = @sqrt(2.0 * std.math.pi);
    return sqrt_2pi * std.math.pow(f64, z + 1.0, z + 0.5) * @exp(-(z + 1.0));
}

pub fn main() !void {
    const allocator = std.heap.page_allocator;
    var args = std.process.args();
    defer args.deinit();

    _ = args.next(); // 跳过程序名
    const filename = args.next() orelse {
        std.debug.print("用法: fft_bin <input_file>\n", .{});
        std.debug.print("功能: 对二进制文件进行FFT分析和GMT0005随机性检测\n", .{});
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

    // 执行GMT0005随机性检测
    std.debug.print("\n开始GMT0005随机性检测...\n", .{});

    var metrics = RandomnessMetrics{
        .monobit_test = monobitTest(buf),
        .frequency_block_test = frequencyBlockTest(buf, 128),
        .runs_test = runsTest(buf),
        .longest_run_test = longestRunTest(buf, 128),
        .spectral_test = try spectralTest(allocator, buf),
    };

    metrics.evaluate();

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
