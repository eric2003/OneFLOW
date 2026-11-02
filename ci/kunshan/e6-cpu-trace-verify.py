#!/usr/bin/env python3
"""Compare OneFLOW main-solver qf1/qf2/invflux CPU traces."""
import argparse
import math
import pathlib
import struct

MAGIC = bytes((79, 70, 84, 82, 67, 48, 49, 0))
HEADER_SIZE = 24

def read_trace(path):
    path = pathlib.Path(path)
    data = path.read_bytes()
    magic, n_faces, n_equations, n_arrays = struct.unpack_from(
        "<8sQII", data, 0)
    if magic != MAGIC or n_arrays != 3:
        raise SystemExit("TRACE_HEADER_FAIL path={}".format(path))
    count = n_faces * n_equations
    offset = HEADER_SIZE
    arrays = []
    for _ in range(3):
        end = offset + count * 8
        if end > len(data):
            raise SystemExit("TRACE_SIZE_FAIL path={}".format(path))
        arrays.append(struct.unpack_from("<{}d".format(count), data, offset))
        offset = end
    if offset != len(data):
        raise SystemExit("TRACE_SIZE_FAIL path={}".format(path))
    return int(n_faces), int(n_equations), arrays

def compare(legacy_path, batch_path):
    legacy = read_trace(legacy_path)
    batch = read_trace(batch_path)
    if legacy[:2] != batch[:2]:
        raise SystemExit(
            "TRACE_SHAPE_FAIL legacy={} batch={}".format(
                legacy[:2], batch[:2]))
    n_faces, n_equations, legacy_arrays = legacy
    _, _, batch_arrays = batch
    overall_absolute = 0.0
    overall_relative = 0.0
    for name, left, right in zip(
            ("qf1", "qf2", "invflux"), legacy_arrays, batch_arrays):
        if not all(math.isfinite(value) for value in left + right):
            raise SystemExit("TRACE_FINITE_FAIL array={}".format(name))
        absolute = max(abs(a - b) for a, b in zip(left, right))
        relative = max(
            abs(a - b) / max(abs(a), abs(b), 1.0e-300)
            for a, b in zip(left, right))
        overall_absolute = max(overall_absolute, absolute)
        overall_relative = max(overall_relative, relative)
        print(
            "TRACE array={} n={} max_absolute={:.17g} "
            "max_relative={:.17g}".format(
                name, len(left), absolute, relative))
    for name, values in (
            ("qf1", legacy_arrays[0]), ("qf2", legacy_arrays[1])):
        density = values[:n_faces]
        pressure = values[4 * n_faces:5 * n_faces]
        if min(density) <= 0.0 or min(pressure) <= 0.0:
            raise SystemExit("TRACE_PHYSICAL_FAIL array={}".format(name))
        print(
            "PHYSICAL array={} finite=true min_density={:.17g} "
            "min_pressure={:.17g}".format(
                name, min(density), min(pressure)))
    print(
        "TRACE_PASS faces={} equations={} max_absolute={:.17g} "
        "max_relative={:.17g}".format(
            n_faces, n_equations, overall_absolute, overall_relative))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("legacy_trace")
    parser.add_argument("batch_trace")
    args = parser.parse_args()
    compare(args.legacy_trace, args.batch_trace)

if __name__ == "__main__":
    main()
