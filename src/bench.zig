//! FFT 基准测试工具 - 基于 zbench 模块
//!
//! 使用 zbench 库对 FFT 各函数进行精确的基准测试，
//! 支持多种数据大小、数据类型和算法变体。
//!
//! 用法:
//!   zig build bench              # 运行所有基准测试
//!   zig build bench -- --size=1024  # 指定 FFT 大小
//!   zig build bench -- --type=f32   # 指定浮点类型

const std = @import("std");
const fft = @import("fft.zig");
const zbench = @import("zbench");

// FFT 库版本（从 build.zig.zon 同步）
const fft_version = "0.1.1";

/// 基准测试配置
pub const BenchConfig = struct {
    /// FFT 大小（必须是 2 的幂）
    size: usize = 1024,
    /// 浮点类型: "f32" 或 "f64" 或 "all"
    float_type: []const u8 = "all",
    /// 是否测试 in-place 版本
    test_inplace: bool = true,
    /// 是否测试 out-of-place 版本
    test_outofplace: bool = true,
    /// 是否测试 DFT（慢速参考）
    test_dft: bool = false,
};

/// 解析命令行参数
fn parseArgs(io: std.Io, args_vec: std.process.Args.Vector, allocator: std.mem.Allocator) !BenchConfig {
    var config = BenchConfig{};
    var args = try std.process.Args.Iterator.initAllocator(
        @as(std.process.Args, .{ .vector = args_vec }),
        allocator,
    );
    defer args.deinit();

    // 跳过程序名
    _ = args.next();

    while (args.next()) |arg| {
        if (std.mem.startsWith(u8, arg, "--size=")) {
            config.size = try std.fmt.parseInt(usize, arg["--size=".len..], 10);
        } else if (std.mem.startsWith(u8, arg, "--type=")) {
            config.float_type = arg["--type=".len..];
        } else if (std.mem.eql(u8, arg, "--no-inplace")) {
            config.test_inplace = false;
        } else if (std.mem.eql(u8, arg, "--no-outofplace")) {
            config.test_outofplace = false;
        } else if (std.mem.eql(u8, arg, "--with-dft")) {
            config.test_dft = true;
        } else if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelp(io) catch {};
            std.process.exit(0);
        }
    }

    return config;
}

fn printHelp(io: std.Io) !void {
    const file: std.Io.File = .stdout();
    var w = file.writerStreaming(io, &.{});
    const writer: *std.Io.Writer = &w.interface;
    try writer.print(
        \\FFT 基准测试工具
        \\
        \\用法: zig build bench -- [选项]
        \\
        \\选项:
        \\  --size=N        FFT 大小（2 的幂，默认: 1024）
        \\  --type=T        浮点类型: f32, f64 或 all（默认: all）
        \\  --no-inplace    跳过 in-place FFT 测试
        \\  --no-outofplace 跳过 out-of-place FFT 测试
        \\  --with-dft      包含 DFT 参考测试（很慢）
        \\  --help, -h      显示此帮助信息
        \\
    , .{});
}

/// 生成测试数据（实数信号）
fn generateRealSignal(comptime T: type, allocator: std.mem.Allocator, size: usize) ![]T {
    const data = try allocator.alloc(T, size);
    for (0..size) |i| {
        const t = @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(size));
        // 混合频率信号: 7Hz + 23Hz 正弦波 + 噪声
        data[i] = @sin(2.0 * std.math.pi * @as(T, 7.0) * t) +
            @as(T, 0.3) * @cos(2.0 * std.math.pi * @as(T, 23.0) * t) +
            @as(T, 0.05) * @sin(2.0 * std.math.pi * @as(T, 101.0) * t);
    }
    return data;
}

/// 生成复数测试数据
fn generateComplexSignal(comptime T: type, allocator: std.mem.Allocator, size: usize) ![]std.math.Complex(T) {
    const data = try allocator.alloc(std.math.Complex(T), size);
    for (0..size) |i| {
        const t = @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(size));
        data[i] = std.math.Complex(T){
            .re = @sin(2.0 * std.math.pi * @as(T, 7.0) * t) + @as(T, 0.3) * @cos(2.0 * std.math.pi * @as(T, 23.0) * t),
            .im = @cos(2.0 * std.math.pi * @as(T, 13.0) * t),
        };
    }
    return data;
}

// ============================================================
// 基准测试函数 - 使用 zbench 的 BenchFunc 签名 (fn(Allocator) void)
// ============================================================

/// f64 in-place FFT 基准测试
fn benchFFTInPlaceF64(allocator: std.mem.Allocator) void {
    const ctx = @as(*struct { size: usize }, @ptrCast(@alignCast(&bench_ctx_f64_inplace)));
    const size = ctx.size;

    // 生成测试数据
    const data = generateComplexSignal(f64, allocator, size) catch @panic("OOM");
    defer allocator.free(data);

    // 复制一份用于测试
    const buf = allocator.alloc(std.math.Complex(f64), size) catch @panic("OOM");
    defer allocator.free(buf);
    @memcpy(buf, data);

    fft.fftInPlace(f64, allocator, buf) catch @panic("FFT failed");
}

/// f64 out-of-place FFT 基准测试
fn benchFFTOutOfPlaceF64(allocator: std.mem.Allocator) void {
    const ctx = @as(*struct { size: usize }, @ptrCast(@alignCast(&bench_ctx_f64_outofplace)));
    const size = ctx.size;

    const input = generateRealSignal(f64, allocator, size) catch @panic("OOM");
    defer allocator.free(input);

    const output = allocator.alloc(std.math.Complex(f64), size) catch @panic("OOM");
    defer allocator.free(output);

    fft.fft(f64, allocator, input, output) catch @panic("FFT failed");
}

/// f32 in-place FFT 基准测试
fn benchFFTInPlaceF32(allocator: std.mem.Allocator) void {
    const ctx = @as(*struct { size: usize }, @ptrCast(@alignCast(&bench_ctx_f32_inplace)));
    const size = ctx.size;

    const data = generateComplexSignal(f32, allocator, size) catch @panic("OOM");
    defer allocator.free(data);

    const buf = allocator.alloc(std.math.Complex(f32), size) catch @panic("OOM");
    defer allocator.free(buf);
    @memcpy(buf, data);

    fft.fftInPlace(f32, allocator, buf) catch @panic("FFT failed");
}

/// f32 out-of-place FFT 基准测试
fn benchFFTOutOfPlaceF32(allocator: std.mem.Allocator) void {
    const ctx = @as(*struct { size: usize }, @ptrCast(@alignCast(&bench_ctx_f32_outofplace)));
    const size = ctx.size;

    const input = generateRealSignal(f32, allocator, size) catch @panic("OOM");
    defer allocator.free(input);

    const output = allocator.alloc(std.math.Complex(f32), size) catch @panic("OOM");
    defer allocator.free(output);

    fft.fft(f32, allocator, input, output) catch @panic("FFT failed");
}

/// DFT f64 基准测试（慢速参考）
fn benchDFTF64(allocator: std.mem.Allocator) void {
    const ctx = @as(*struct { size: usize }, @ptrCast(@alignCast(&bench_ctx_dft)));
    const size = ctx.size;

    const input = generateComplexSignal(f64, allocator, size) catch @panic("OOM");
    defer allocator.free(input);

    const output = allocator.alloc(std.math.Complex(f64), size) catch @panic("OOM");
    defer allocator.free(output);

    fft.dft(f64, input, output);
}

// 全局上下文变量（用于传递参数给 BenchFunc）
var bench_ctx_f64_inplace: struct { size: usize } = .{ .size = 1024 };
var bench_ctx_f64_outofplace: struct { size: usize } = .{ .size = 1024 };
var bench_ctx_f32_inplace: struct { size: usize } = .{ .size = 1024 };
var bench_ctx_f32_outofplace: struct { size: usize } = .{ .size = 1024 };
var bench_ctx_dft: struct { size: usize } = .{ .size = 1024 };

/// 注册所有基准测试到 zbench
fn registerBenchmarks(bench: *zbench.Benchmark, config: BenchConfig) !void {
    const allocator = bench.allocator;
    const sizes = [_]usize{ 64, 256, 1024, 4096, 16384 };

    // 根据配置决定测试哪些类型
    const test_f64 = std.mem.eql(u8, config.float_type, "f64") or std.mem.eql(u8, config.float_type, "all");
    const test_f32 = std.mem.eql(u8, config.float_type, "f32") or std.mem.eql(u8, config.float_type, "all");

    // 多尺寸基准测试
    for (sizes) |size| {
        if (test_f64) {
            if (config.test_inplace) {
                bench_ctx_f64_inplace.size = size;
                const name = try std.fmt.allocPrint(allocator, "I f64 {d}", .{size});
                try bench.add(name, benchFFTInPlaceF64, .{});
            }
            if (config.test_outofplace) {
                bench_ctx_f64_outofplace.size = size;
                const name = try std.fmt.allocPrint(allocator, "O f64 {d}", .{size});
                try bench.add(name, benchFFTOutOfPlaceF64, .{});
            }
        }
        if (test_f32) {
            if (config.test_inplace) {
                bench_ctx_f32_inplace.size = size;
                const name = try std.fmt.allocPrint(allocator, "I f32 {d}", .{size});
                try bench.add(name, benchFFTInPlaceF32, .{});
            }
            if (config.test_outofplace) {
                bench_ctx_f32_outofplace.size = size;
                const name = try std.fmt.allocPrint(allocator, "O f32 {d}", .{size});
                try bench.add(name, benchFFTOutOfPlaceF32, .{});
            }
        }
    }

    // 单尺寸详细测试
    if (test_f64) {
        if (config.test_inplace) {
            bench_ctx_f64_inplace.size = config.size;
            const name = try std.fmt.allocPrint(allocator, "I f64 {d}", .{config.size});
            try bench.add(name, benchFFTInPlaceF64, .{ .iterations = 100, .max_iterations = 500 });
        }
        if (config.test_outofplace) {
            bench_ctx_f64_outofplace.size = config.size;
            const name = try std.fmt.allocPrint(allocator, "O f64 {d}", .{config.size});
            try bench.add(name, benchFFTOutOfPlaceF64, .{ .iterations = 100, .max_iterations = 500 });
        }
        if (config.test_dft) {
            bench_ctx_dft.size = config.size;
            const name = try std.fmt.allocPrint(allocator, "DFT {d}", .{config.size});
            try bench.add(name, benchDFTF64, .{ .iterations = 3, .max_iterations = 10 });
        }
    }

    if (test_f32) {
        if (config.test_inplace) {
            bench_ctx_f32_inplace.size = config.size;
            const name = try std.fmt.allocPrint(allocator, "I f32 {d}", .{config.size});
            try bench.add(name, benchFFTInPlaceF32, .{ .iterations = 100, .max_iterations = 500 });
        }
        if (config.test_outofplace) {
            bench_ctx_f32_outofplace.size = config.size;
            const name = try std.fmt.allocPrint(allocator, "O f32 {d}", .{config.size});
            try bench.add(name, benchFFTOutOfPlaceF32, .{ .iterations = 100, .max_iterations = 500 });
        }
    }
}

// ============================================================
// 主入口
// ============================================================

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const stdout: std.Io.File = .stdout();
    var stdout_w = stdout.writerStreaming(io, &.{});
    const writer: *std.Io.Writer = &stdout_w.interface;

    const config = try parseArgs(io, init.minimal.args.vector, init.gpa);

    // 获取当前日期时间
    const now_ts = std.Io.Timestamp.now(io, .real);
    const epoch_seconds: i64 = @intCast(@divTrunc(now_ts.nanoseconds, @as(i96, @intCast(std.time.ns_per_s))));
    const epoch_day = std.time.epoch.EpochSeconds{ .secs = @intCast(epoch_seconds) };
    const year_and_day = epoch_day.getEpochDay().calculateYearDay();
    const month_and_day = year_and_day.calculateMonthDay();
    const day_seconds = epoch_day.getDaySeconds();

    // 获取系统信息（CPU、OS等）
    const sysinfo = zbench.getSystemInfo() catch |err| blk: {
        writer.print("[警告] 无法获取系统信息: {}\n", .{err}) catch {};
        break :blk null;
    };

    // 获取编译优化模式
    const optimize_mode = @import("builtin").mode;
    const mode_str = switch (optimize_mode) {
        .Debug => "Debug",
        .ReleaseSafe => "ReleaseSafe",
        .ReleaseFast => "ReleaseFast",
        .ReleaseSmall => "ReleaseSmall",
    };

    writer.print(
        \\╔══════════════════════════════════════════╗
        \\║      FFT 基准测试工具 (zbench)           ║
        \\╚══════════════════════════════════════════╝

        \\ I: In-Place
        \\ O: Out-of-Place
        \\
        \\测试时间: {d:0>4}-{d:0>2}-{d:0>2} {d:0>2}:{d:0>2}:{d:0>2}
        \\优化模式: {s}
        \\FFT 版本: {s}
        \\
        \\配置:
        \\  FFT 大小: {d}
        \\  浮点类型: {s}
        \\  In-Place: {s}
        \\  Out-of-Place: {s}
        \\  DFT 参考: {s}
        \\
    , .{
        year_and_day.year,
        month_and_day.month.numeric(),
        month_and_day.day_index + 1,
        day_seconds.getHoursIntoDay(),
        day_seconds.getMinutesIntoHour(),
        day_seconds.getSecondsIntoMinute(),
        mode_str,
        fft_version,
        config.size,
        config.float_type,
        if (config.test_inplace) "yes" else "no",
        if (config.test_outofplace) "yes" else "no",
        if (config.test_dft) "yes" else "no",
    }) catch {};

    // 打印系统信息
    if (sysinfo) |info| {
        writer.print(
            \\
            \\系统信息:
            \\  操作系统: {s}
            \\  CPU:      {s}
            \\  核心数:   {d}
            \\  内存:     {Bi:.3}
            \\
        , .{
            info.platform,
            info.cpu,
            info.cpu_cores,
            info.memory_total,
        }) catch {};
    }

    // 初始化 zbench
    var bench = zbench.Benchmark.init(init.gpa, .{});
    defer bench.deinit();

    // 注册所有基准测试
    try registerBenchmarks(&bench, config);

    writer.print("\n", .{}) catch {};
    writer.print("共 {d} 个基准测试\n\n", .{bench.benchmarks.items.len}) catch {};

    // 运行基准测试
    try bench.run(io, stdout);

    // 释放 benchmark 名称（zbench 只保存引用，不拥有所有权）
    for (bench.benchmarks.items) |item| {
        bench.allocator.free(item.name);
    }

    writer.print("\n基准测试完成！\n", .{}) catch {};
}
