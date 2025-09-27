//! 基础类型与SIMD支持
const std = @import("std");

pub const Complex = std.math.Complex(f64);
pub const VectorF64 = @Vector(4, f64);
pub const VectorF64x8 = @Vector(8, f64);

pub const VectorComplex = struct {
    re: VectorF64,
    im: VectorF64,
};
