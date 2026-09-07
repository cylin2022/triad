"""
resource_tuner.py
Hardware-Aware Auto-Tuning Engine for TRIAD Platform.
Dynamically profiles host CPU cores, memory limits, RAM disk (/dev/shm), and system load
to maximize parallel throughput and eliminate I/O latency.
"""

import os
import shutil
import platform
import multiprocessing
from pathlib import Path
from typing import Dict, Any, Tuple
import psutil

# ─── RAM Disk (/dev/shm) Optimization ─────────────────────────────────────────

def get_ram_disk_tmpdir() -> Path:
    """
    Returns high-speed in-memory RAM disk directory (/dev/shm/netmhc_tmp).
    Eliminates physical disk I/O bottlenecks for temporary peptide and FASTA files.
    """
    shm_path = Path("/dev/shm/netmhc_tmp")
    try:
        if Path("/dev/shm").exists() and os.access("/dev/shm", os.W_OK):
            shm_path.mkdir(parents=True, exist_ok=True)
            # Test write
            test_file = shm_path / ".tuner_test"
            test_file.write_text("ok")
            test_file.unlink(missing_ok=True)
            return shm_path
    except Exception:
        pass
    
    # Fallback to local project tmp
    fallback = Path(__file__).parent.parent / "netMHCpan-4.2" / "tmp"
    fallback.mkdir(parents=True, exist_ok=True)
    return fallback


def cleanup_ram_disk():
    """Clean up stale temp files in RAM disk."""
    try:
        shm_path = Path("/dev/shm/netmhc_tmp")
        if shm_path.exists():
            for f in shm_path.glob("netMHCpan_*"):
                if f.is_dir():
                    shutil.rmtree(f, ignore_errors=True)
                else:
                    f.unlink(missing_ok=True)
            for f in shm_path.glob("*.pep"):
                f.unlink(missing_ok=True)
            for f in shm_path.glob("*.fsa"):
                f.unlink(missing_ok=True)
    except Exception:
        pass


# ─── Hardware Profiler ────────────────────────────────────────────────────────

def get_cpu_model_name() -> str:
    """Extract human-readable CPU model string."""
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if "model name" in line:
                    return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return platform.processor() or "Multi-Core x86_64 CPU"


def get_system_specs() -> Dict[str, Any]:
    """
    Returns comprehensive hardware specs and real-time utilization metrics.
    """
    logical_cpus = psutil.cpu_count(logical=True) or multiprocessing.cpu_count()
    physical_cpus = psutil.cpu_count(logical=False) or (logical_cpus // 2)
    load_avg = os.getloadavg() if hasattr(os, "getloadavg") else (0.0, 0.0, 0.0)
    cpu_percent = psutil.cpu_percent(interval=None)

    # Memory
    mem = psutil.virtual_memory()
    total_ram_gb = round(mem.total / (1024 ** 3), 1)
    avail_ram_gb = round(mem.available / (1024 ** 3), 1)
    used_ram_gb  = round(mem.used / (1024 ** 3), 1)
    ram_percent  = mem.percent

    # RAM Disk (/dev/shm)
    ram_disk_active = False
    shm_total_gb = 0.0
    shm_avail_gb = 0.0
    if Path("/dev/shm").exists():
        try:
            shm_usage = shutil.disk_usage("/dev/shm")
            shm_total_gb = round(shm_usage.total / (1024 ** 3), 1)
            shm_avail_gb = round(shm_usage.free / (1024 ** 3), 1)
            ram_disk_active = True
        except Exception:
            pass

    # Disk (Root /)
    root_usage = shutil.disk_usage("/")
    root_total_tb = round(root_usage.total / (1024 ** 4), 2)
    root_free_tb  = round(root_usage.free / (1024 ** 4), 2)
    root_used_pct = round((root_usage.used / root_usage.total) * 100, 1)

    # Dynamic parallel worker pool capacity
    optimal_workers = max(2, min(logical_cpus - 2, 32))

    return {
        "cpu": {
            "model": get_cpu_model_name(),
            "logical_cores": logical_cpus,
            "physical_cores": physical_cpus,
            "load_1m": round(load_avg[0], 2),
            "load_5m": round(load_avg[1], 2),
            "load_15m": round(load_avg[2], 2),
            "cpu_percent": cpu_percent,
            "optimal_workers": optimal_workers,
        },
        "memory": {
            "total_gb": total_ram_gb,
            "available_gb": avail_ram_gb,
            "used_gb": used_ram_gb,
            "percent": ram_percent,
        },
        "ram_disk": {
            "active": ram_disk_active,
            "path": "/dev/shm/netmhc_tmp",
            "total_gb": shm_total_gb,
            "available_gb": shm_avail_gb,
            "description": "Zero-latency In-Memory Temp Storage",
        },
        "disk": {
            "total_tb": root_total_tb,
            "free_tb": root_free_tb,
            "used_percent": root_used_pct,
        },
    }


# ─── Dynamic Batch Configuration ──────────────────────────────────────────────

def get_optimal_batch_config(n_items: int, n_alleles: int = 1) -> Dict[str, Any]:
    """
    Computes optimal CPU worker allocation and chunk size based on dataset scale
    and real-time host hardware capabilities.
    """
    specs = get_system_specs()
    logical_cpus = specs["cpu"]["logical_cores"]
    current_load = specs["cpu"]["load_1m"]

    # Compute worker scaling: if system load is high, throttle slightly; if idle, scale up
    load_ratio = current_load / max(1, logical_cpus)
    if load_ratio > 0.85:
        max_workers = max(2, logical_cpus // 4)
    elif load_ratio > 0.60:
        max_workers = max(4, logical_cpus // 3)
    else:
        max_workers = max(4, min(logical_cpus - 2, 32))

    # Adaptive chunk size
    if n_items <= 50:
        chunk_size = max(5, n_items // 4 if n_items >= 10 else n_items)
    elif n_items <= 500:
        chunk_size = max(20, n_items // max_workers)
    elif n_items <= 5000:
        chunk_size = 150
    else:
        chunk_size = 400

    n_chunks = max(1, (n_items + chunk_size - 1) // chunk_size)
    allocated_workers = min(n_chunks, max_workers)

    return {
        "allocated_workers": allocated_workers,
        "chunk_size": chunk_size,
        "n_chunks": n_chunks,
        "ram_disk_active": specs["ram_disk"]["active"],
        "summary": (
            f"⚡ Hardware Auto-Tuning: {allocated_workers} CPU Workers across {logical_cpus} cores "
            f"| Chunk Size: {chunk_size} | RAM Disk (/dev/shm): Active"
        ),
    }


if __name__ == "__main__":
    import pprint
    print("=== Host Hardware Specs ===")
    pprint.pprint(get_system_specs())
    print("\n=== Sample Batch Config for 1,000 peptides ===")
    pprint.pprint(get_optimal_batch_config(1000, 2))
