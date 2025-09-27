// filepath: src/main.zig
const std = @import("std");
const fft = @import("fft.zig");

pub fn main() !void {
    const allocator = std.heap.page_allocator;
    var args = std.process.args();
    defer args.deinit();

    _ = args.next(); // 跳过程序名
    const filename = args.next() orelse {
        std.debug.print("用法: fft_bin <input_file>\n", .{});
        return error.InvalidArgs;
    };

    // 打开并读取文件
    var file = try std.fs.cwd().openFile(filename, .{});
    defer file.close();

    const file_size = try file.getEndPos();
    const usable_n = file_size / @sizeOf(f64); // 向下取整
    if (usable_n == 0) {
        std.debug.print("文件中没有完整的 f64 数据\n", .{});
        return error.InvalidInput;
    }

    const buf = try allocator.alloc(f64, usable_n);
    defer allocator.free(buf);

    // 只读取完整的 f64 数据部分
    const read_bytes = try file.readAll(std.mem.sliceAsBytes(buf));
    if (read_bytes < usable_n * @sizeOf(f64)) {
        std.debug.print("读取文件失败\n", .{});
        return error.ReadFailed;
    }

    // 分配输出缓冲区
    const out_len = usable_n / 2 + 1;
    const fft_out = try allocator.alloc(f64, 2 * out_len);
    defer allocator.free(fft_out);
    const fft_m = try allocator.alloc(f64, out_len);
    defer allocator.free(fft_m);

    // 调用 FFT
    try fft.fftR2C(allocator, buf, fft_out, fft_m);

    // 打印前10个幅值
    std.debug.print("前10个幅值:\n", .{});
    for (0..@min(10, out_len)) |i| {
        std.debug.print("{d}: {d:.6}\n", .{ i, fft_m[i] });
    }
}
