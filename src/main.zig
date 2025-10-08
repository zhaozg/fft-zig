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

    // 分配输出缓冲区
    const out_len = nbits / 2 + 1;
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
