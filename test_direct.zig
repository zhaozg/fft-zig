const std = @import("std");
const fft = @import("src/fft.zig");
const Complex = fft.Complex;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();
    
    const size: usize = 1024;
    const freq: usize = 10;
    
    const data = try allocator.alloc(Complex, size);
    defer allocator.free(data);
    
    for (0..size) |i| {
        const t = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(size));
        data[i] = Complex{ 
            .re = @sin(2.0 * std.math.pi * @as(f64, @floatFromInt(freq)) * t),
            .im = 0.0
        };
    }
    
    try fft.fftInPlace(allocator, data);
    
    const mag = @sqrt(data[freq].re * data[freq].re + data[freq].im * data[freq].im);
    std.debug.print("Size: {d}, Magnitude at bin {d}: {d:.4}\n", .{size, freq, mag});
    std.debug.print("Expected: {d:.4}\n", .{@as(f64, @floatFromInt(size)) / 2.0});
    std.debug.print("Ratio: {d:.6}\n", .{mag / (@as(f64, @floatFromInt(size)) / 2.0)});
}
