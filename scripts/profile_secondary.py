"""Bounded CPU sizing probe; synthetic token IDs are NEVER biological training data.

Runs each configuration in a fresh, below-normal-priority Windows process. This
is engineering profiling, not a hyperparameter trial or a biological experiment.
The JSON plan is persisted before the first probe; no checkpoint is produced.
"""
from __future__ import annotations

import argparse
import ctypes
from ctypes import wintypes
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
GIB = 1024**3
RESERVE_BYTES = GIB
PLAN = [(48, 2, 4, 4), (96, 3, 4, 4), (192, 4, 6, 4),
        (192, 4, 6, 8), (256, 6, 8, 4), (384, 6, 8, 4)]


def physical_memory():
    class MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [("dwLength", wintypes.DWORD), ("dwMemoryLoad", wintypes.DWORD)] + [
            (field, ctypes.c_ulonglong) for field in
            ("ullTotalPhys", "ullAvailPhys", "ullTotalPageFile", "ullAvailPageFile",
             "ullTotalVirtual", "ullAvailVirtual", "ullAvailExtendedVirtual")]
    info = MEMORYSTATUSEX()
    info.dwLength = ctypes.sizeof(info)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(info)):
        raise ctypes.WinError()
    return {"total_physical_bytes": info.ullTotalPhys,
            "available_physical_bytes": info.ullAvailPhys}


def process_memory():
    class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [
            (field, ctypes.c_size_t) for field in
            ("PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage",
             "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage", "QuotaNonPagedPoolUsage",
             "PagefileUsage", "PeakPagefileUsage")]
    kernel = ctypes.windll.kernel32
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    read_memory = ctypes.windll.psapi.GetProcessMemoryInfo
    read_memory.argtypes = [wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD]
    read_memory.restype = wintypes.BOOL
    info = PROCESS_MEMORY_COUNTERS()
    info.cb = ctypes.sizeof(info)
    if not read_memory(
            kernel.GetCurrentProcess(), ctypes.byref(info), info.cb):
        raise ctypes.WinError()
    return {"rss_bytes": info.WorkingSetSize,
            "peak_rss_bytes": info.PeakWorkingSetSize,
            "private_commit_bytes": info.PagefileUsage}


def estimates(width, layers, heads, batch=16, tokens=100, vocab_size=64):
    parameters = vocab_size * width + layers * (12 * width * width + 13 * width) + 2 * width
    # Parameter, gradient and two FP32 Adam moments, with two extra parameter
    # copies as headroom; activation estimate explicitly includes quadratic
    # attention. This is a conservative planning heuristic, not a proof/bound.
    parameter_bytes = parameters * 24
    activation_bytes = 4 * batch * layers * (
        18 * tokens * width + 4 * heads * tokens * tokens)
    return {"estimated_parameters": parameters,
            "estimated_training_allocation_bytes": parameter_bytes + activation_bytes,
            "estimated_process_bytes": parameter_bytes + activation_bytes + 384 * 1024**2}


def spec_estimates(spec):
    return estimates(spec["d_model"], spec["n_layers"], spec["n_heads"],
                     spec["batch_size"], spec["tokens"], spec.get("vocab_size", 64))


def child(spec):
    import torch
    sys.path.insert(0, str(ROOT / "src"))
    from microprotein_lm.model import Decoder, ModelConfig, summed_loss

    torch.set_num_threads(spec["threads"])
    torch.set_num_interop_threads(1)
    torch.manual_seed(20260925)
    mem = physical_memory()
    estimate = spec_estimates(spec)
    result = {"config": spec, "environment": {
        "python": platform.python_version(), "torch": torch.__version__,
        "cuda_available": torch.cuda.is_available(), "xpu_available": torch.xpu.is_available(),
        "platform": platform.platform(), "cpu_logical_count": os.cpu_count(),
        "dtype": "float32", "device": "cpu", "interop_threads": 1},
        "memory_after_import": mem, "rss_after_import": process_memory(), **estimate}
    if mem["available_physical_bytes"] < RESERVE_BYTES + estimate["estimated_training_allocation_bytes"]:
        return {**result, "status": "skipped_memory_gate_after_import"}
    vocab_size = spec.get("vocab_size", 64)
    config = ModelConfig(vocab_size=vocab_size, max_tokens=spec["tokens"], token_width=3 if vocab_size == 64 else 1,
        d_model=spec["d_model"], n_layers=spec["n_layers"],
        n_heads=spec["n_heads"], dropout=0.1)
    model = Decoder(config)
    model.train()
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.001, weight_decay=0.01, foreach=False)
    ids = torch.randint(vocab_size, (spec["batch_size"], spec["tokens"] + 1))
    x, targets = ids[:, :-1], ids[:, 1:].clone()
    if vocab_size == 4:
        targets[:, :2] = -100
    target_count = (targets != -100).sum().item()
    durations = []
    minimum_available = physical_memory()["available_physical_bytes"]
    for step in range(7):
        available = physical_memory()["available_physical_bytes"]
        minimum_available = min(minimum_available, available)
        if available < RESERVE_BYTES:
            return {**result, "status": "stopped_desktop_memory_reserve",
                    "completed_steps": step, "minimum_available_physical_bytes": minimum_available,
                    "process_memory": process_memory()}
        started = time.perf_counter()
        optimizer.zero_grad(set_to_none=True)
        logits = model(x)
        loss = summed_loss(logits, targets) / target_count
        loss.backward()
        norm = torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        elapsed = time.perf_counter() - started
        if not torch.isfinite(loss) or not torch.isfinite(norm):
            raise RuntimeError("Nonfinite synthetic sizing probe")
        if step >= 2:
            durations.append(elapsed)
        minimum_available = min(minimum_available, physical_memory()["available_physical_bytes"])
    median = statistics.median(durations)
    return {**result, "status": "completed", "parameters": sum(p.numel() for p in model.parameters()),
            "warmup_steps": 2, "measured_steps": 5, "step_seconds": durations,
            "median_step_seconds": median, "min_step_seconds": min(durations),
            "max_step_seconds": max(durations),
            "estimated_1200_update_seconds_excluding_diagnostics": median * 1200,
            "estimated_2000_update_seconds_excluding_diagnostics": median * 2000,
            "minimum_available_physical_bytes": minimum_available,
            "process_memory": process_memory()}


def run_spec(spec):
    mem = physical_memory()
    estimate = spec_estimates(spec)
    if mem["available_physical_bytes"] < RESERVE_BYTES + estimate["estimated_process_bytes"]:
        return {"config": spec, "status": "skipped_memory_gate_before_import",
                "host_memory": mem, **estimate}
    print(f"Sizing {spec}; available RAM {mem['available_physical_bytes']/GIB:.2f} GiB", flush=True)
    started = time.perf_counter()
    try:
        proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--child", json.dumps(spec)],
            cwd=ROOT, capture_output=True, text=True, timeout=90,
            creationflags=subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        if proc.returncode:
            result = {"config": spec, "status": "process_error", "stderr": proc.stderr[-4000:]}
        else:
            result = json.loads(proc.stdout)
    except subprocess.TimeoutExpired:
        result = {"config": spec, "status": "process_timeout"}
    result["process_wall_seconds"] = time.perf_counter() - started
    return result


def extend_matrix(output):
    report = json.loads(output.read_text(encoding="utf-8"))
    if "matrix_followup" in report:
        raise SystemExit("Refusing to overwrite existing matrix followup")
    plan = [{"d_model": d, "n_layers": n, "n_heads": h, "threads": 4,
             "batch_size": 16, "tokens": tokens, "vocab_size": vocab}
            for d, n, h, tokens, vocab in [(384, 6, 8, 206, 4), (384, 6, 8, 302, 4),
                (384, 6, 8, 68, 64), (192, 4, 6, 206, 4), (192, 4, 6, 302, 4), (192, 4, 6, 68, 64)]]
    followup = {"created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "planned_before_execution", "plan": plan,
        "purpose": "Profile matched base contexts and realistic ATP8 codon length; retain 192x4 fallback if 384x6 base contexts fail conservative memory gate",
        "protocol": report["protocol"], "host_memory_before": physical_memory(),
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "model_sha256": hashlib.sha256((ROOT / "src/microprotein_lm/model.py").read_bytes()).hexdigest(),
        "results": []}
    report["matrix_followup"] = followup
    def save():
        output.write_bytes((json.dumps(report, indent=2) + "\n").encode("utf-8"))
    save()
    for spec in plan:
        result = run_spec(spec)
        followup["results"].append(result)
        followup["status"] = "running"
        save()
        print(json.dumps({k: result[k] for k in ("config", "status", "median_step_seconds", "process_memory") if k in result}), flush=True)
    all_results = report["results"] + followup["results"]
    def find(width, layers, vocab, tokens):
        return next((r for r in all_results if r["status"] == "completed"
            and r["config"]["d_model"] == width and r["config"]["n_layers"] == layers
            and r["config"].get("vocab_size", 64) == vocab and r["config"]["tokens"] == tokens
            and r["config"]["threads"] == 4), None)
    alternatives = []
    small = find(48, 2, 64, 100)
    for width, layers in [(384, 6), (192, 4)]:
        base = find(width, layers, 4, 206)
        codon = find(width, layers, 64, 100)
        if base and codon and small:
            # 8 arms x 3 seeds + 1 full-pool run. One arm is base, one is the
            # original small codon decoder, six arms plus full-pool are large codon.
            seconds_per_matrix_update = (3 * base["median_step_seconds"] +
                19 * codon["median_step_seconds"] + 3 * small["median_step_seconds"])
            budget_limit = int(4 * 3600 / (1.25 * seconds_per_matrix_update))
            recommendation = max(100, int(budget_limit * 0.8 // 100) * 100)
            alternatives.append({"d_model": width, "n_layers": layers,
                "batch_size": 16, "threads": 4, "matrix_runs": 25,
                "seconds_per_one_update_in_all_runs": seconds_per_matrix_update,
                "maximum_updates_from_short_timing_with_25_percent_overhead": budget_limit,
                "recommended_updates_per_run": recommendation,
                "estimated_matrix_seconds_with_25_percent_overhead": recommendation * seconds_per_matrix_update * 1.25,
                "assumptions": "3 base runs at T206, 19 large codon runs charged at worst-case T100, 3 small codon runs at T100; effective batch16; sequential execution. Recommendation holds back another20% for sustained slowdown. Full-pool diagnostics must fit within overhead or be separately budgeted; this is an estimate, not a guaranteed wall-time limit."})
    followup["matrix_budget_alternatives"] = alternatives
    followup["status"] = "completed"
    followup["completed_utc"] = datetime.now(timezone.utc).isoformat()
    followup["host_memory_after"] = physical_memory()
    save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child", help=argparse.SUPPRESS)
    parser.add_argument("--output", default="reports/secondary/hardware-benchmark.json")
    parser.add_argument("--extend-matrix", action="store_true", help="Append prespecified matched-base followup without replacing initial results")
    args = parser.parse_args()
    if platform.system() != "Windows":
        raise SystemExit("This probe uses Windows native memory measurements.")
    if args.child:
        print(json.dumps(child(json.loads(args.child))))
        return
    output = ROOT / args.output
    if args.extend_matrix:
        extend_matrix(output)
        return
    if output.exists():
        raise SystemExit(f"Refusing to overwrite prior profile: {output}; choose a new --output")
    output.parent.mkdir(parents=True, exist_ok=True)
    plan = [{"d_model": d, "n_layers": n, "n_heads": h, "threads": threads,
             "batch_size": 16, "tokens": 100} for d, n, h, threads in PLAN]
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "Synthetic hardware profiling only; no biological examples, scientific training or checkpoints",
        "status": "planned_before_execution", "plan": plan,
        "protocol": {"warmup_steps": 2, "measured_steps": 5,
            "isolated_process_per_config": True, "process_priority": "below_normal",
            "physical_memory_reserve_bytes": RESERVE_BYTES,
            "per_process_timeout_seconds": 90, "optimizer": "AdamW, foreach=False",
            "gradient_clip_norm": 1.0, "seed": 20260925,
            "selection_rule": "Largest completed candidate with >=1 GiB physical RAM available and median step <2 seconds",
            "limitations": "Short synthetic profiling; not sustained thermals, absolute hardware ceiling, data loading, checkpoint or full-corpus diagnostic timing. No concurrent training jobs."},
        "host_memory_before": physical_memory(),
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "model_sha256": hashlib.sha256((ROOT / "src/microprotein_lm/model.py").read_bytes()).hexdigest(),
        "gpu_assessment": {
            "installed_torch": "CPU-only build; directly checked xpu.is_available() is False",
            "observed_hardware_from_parent_inventory": {"cpu": "Intel Core i7-12700H",
                "physical_cores": 14, "logical_processors": 20,
                "gpu": [{"name": "Intel Arc A370M", "driver": "31.0.101.4824"},
                        {"name": "Intel Iris Xe", "driver": "31.0.101.4502"}]},
            "finding": "PyTorch 2.14 documents Arc A-Series support on Windows 11. Intel 2.14 prerequisites specify Windows driver 32.0.101.8801 or newer. Installed CPU build and older driver do not establish a usable XPU environment. Iris Xe is not listed among validated client GPUs on these pages. No software or driver change attempted; no Arc training benchmark performed.",
            "source_access_date_utc": datetime.now(timezone.utc).date().isoformat(),
            "sources": ["https://docs.pytorch.org/docs/2.14/notes/get_start_xpu.html",
                        "https://www.intel.com/content/www/us/en/developer/articles/tool/pytorch-prerequisites-for-intel-gpu/2-14.html"]},
        "results": []}
    def save():
        output.write_bytes((json.dumps(report, indent=2) + "\n").encode("utf-8"))
    save()
    for spec in plan:
        result = run_spec(spec)
        report["results"].append(result)
        report["status"] = "running"
        save()
        print(json.dumps({k: result[k] for k in ("config", "status", "median_step_seconds", "process_memory") if k in result}), flush=True)
    feasible = [r for r in report["results"] if r["status"] == "completed"
                and r["minimum_available_physical_bytes"] >= RESERVE_BYTES
                and r["median_step_seconds"] < 2]
    report["recommended_largest_profiled_config"] = max(feasible, key=lambda r: (
        r["parameters"], -r["median_step_seconds"]))["config"] if feasible else None
    report["status"] = "completed"
    report["completed_utc"] = datetime.now(timezone.utc).isoformat()
    report["host_memory_after"] = physical_memory()
    save()


if __name__ == "__main__":
    main()
