//! 基础类型与SIMD支持
const std = @import("std");

/// 泛型复数类型
pub fn Complex(comptime T: type) type {
    return std.math.Complex(T);
}

/// 泛型向量类型
pub fn Vector(comptime T: type, comptime len: usize) type {
    return @Vector(len, T);
}

/// 4元素向量（用于SIMD优化）
pub fn Vector4(comptime T: type) type {
    return Vector(T, 4);
}

/// 8元素向量（用于SIMD优化）
pub fn Vector8(comptime T: type) type {
    return Vector(T, 8);
}

/// 复数向量结构体（用于SIMD优化）
pub fn VectorComplex(comptime T: type) type {
    return struct {
        re: Vector4(T),
        im: Vector4(T),
    };
}

// 向后兼容的类型别名
pub const ComplexF64 = Complex(f64);
pub const VectorF64 = Vector4(f64);
pub const VectorF64x8 = Vector8(f64);
pub const VectorComplexF64 = VectorComplex(f64);
