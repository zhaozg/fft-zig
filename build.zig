const std = @import("std");

fn getRunArgs(b: *std.Build) ?[]const []const u8 {
    if (@hasField(std.Build, "args")) {
        return b.args;
    } else {
        return b.option([]const []const u8, "args", "Arguments to pass to the executable");
    }
}

pub fn module(b: *std.Build) *std.Build.Module {
    return b.createModule(.{
        .root_source_file = b.path("src/fft.zig"),
    });
}

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const run_args = getRunArgs(b);

    const lib_mod = b.createModule(.{
        .root_source_file = b.path("src/fft.zig"),
        .target = target,
        .optimize = optimize,
    });

    const exe_mod = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });
    exe_mod.addImport("fft", lib_mod);

    const lib = b.addLibrary(.{
        .linkage = .static,
        .name = "fft",
        .root_module = lib_mod,
    });
    b.installArtifact(lib);

    const exe = b.addExecutable(.{
        .name = "fft",
        .root_module = exe_mod,
    });
    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (run_args) |args| run_cmd.addArgs(args);
    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    // ============================================================
    // 基准测试 (zbench) – 直接使用源码，不依赖其 build.zig 的模块
    // ============================================================

    const zbench_dep = b.dependency("zbench", .{
        .target = target,
        .optimize = optimize,
        .with_examples = false, // 跳过示例（避免 addModule 冲突）
    });

    // 直接从依赖包中取 zbench 的根源文件
    const zbench_mod = b.createModule(.{
        .root_source_file = zbench_dep.path("src/zbench.zig"),
        .target = target,
        .optimize = optimize,
    });

    const bench_mod = b.createModule(.{
        .root_source_file = b.path("src/bench.zig"),
        .target = target,
        .optimize = optimize,
    });
    bench_mod.addImport("fft", lib_mod);
    bench_mod.addImport("zbench", zbench_mod);

    const bench_exe = b.addExecutable(.{
        .name = "bench",
        .root_module = bench_mod,
    });
    b.installArtifact(bench_exe);

    const bench_cmd = b.addRunArtifact(bench_exe);
    bench_cmd.step.dependOn(b.getInstallStep());
    if (run_args) |args| bench_cmd.addArgs(args);
    const bench_step = b.step("bench", "Run FFT benchmarks using zbench");
    bench_step.dependOn(&bench_cmd.step);

    // ============================================================
    // 测试
    // ============================================================

    const tst_mod = b.createModule(.{
        .root_source_file = b.path("src/test.zig"),
        .target = target,
        .optimize = optimize,
    });
    tst_mod.addImport("fft", lib_mod);

    const lib_unit_tests = b.addTest(.{ .root_module = tst_mod });
    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);
    const exe_unit_tests = b.addTest(.{ .root_module = tst_mod });
    const run_exe_unit_tests = b.addRunArtifact(exe_unit_tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);
    test_step.dependOn(&run_exe_unit_tests.step);
}
