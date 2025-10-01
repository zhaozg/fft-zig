const std = @import("std");
const fft = @import("src/fft.zig");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();
    
    const size: usize = 1048576;
    const freq: usize = 5;
    
    var input = try allocator.alloc(f64, size);
    defer allocator.free(input);
    
    for (0..size) |i| {
        const t = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(size));
        input[i] = @sin(2.0 * std.math.pi * @as(f64, @floatFromInt(freq)) * t);
    }
    
    const output = try allocator.alloc(f64, size + 2);
    defer allocator.free(output);
    
    const magnitude = try allocator.alloc(f64, size / 2 + 1);
    defer allocator.free(magnitude);
    
    try fft.fftR2C(allocator, input, output, magnitude);
    
    std.debug.print("Magnitude at bin {d}: {d:.4}\n", .{freq, magnitude[freq]});
    std.debug.print("Expected: {d:.4}\n", .{@as(f64, @floatFromInt(size)) / 2.0});
    std.debug.print("Ratio: {d:.6}\n", .{magnitude[freq] / (@as(f64, @floatFromInt(size)) / 2.0)});
}
