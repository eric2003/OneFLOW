#!/usr/bin/env python3
"""Streaming verifier for legacy CPU, CPU batch, and HIP batch stage traces."""
import argparse
import math
import pathlib
import struct


MAGIC = bytes((79, 70, 83, 84, 71, 48, 49, 0))
HEADER = struct.Struct("<8sIIiiIIQ")
KIND_NAMES = {1: "face", 2: "inviscid_residual", 3: "state", 4: "face_metadata"}
ARRAY_NAMES = {
    1: ("qf1", "qf2", "invflux"),
    2: ("inviscid_residual",),
    3: ("state",),
    4: ("left_cell", "right_cell", "boundary_mask", "normal_x", "normal_y",
        "normal_z", "mesh_velocity_normal", "face_area"),
}


def fail(message):
    raise SystemExit(message)


def read_exact(stream, size, path):
    payload = stream.read(size)
    if len(payload) != size:
        fail("STAGE_TRACE_SIZE_FAIL path={}".format(path))
    return payload


def iter_trace(path):
    path = pathlib.Path(path)
    count_records = 0
    with path.open("rb") as stream:
        while True:
            payload = stream.read(HEADER.size)
            if not payload:
                break
            if len(payload) != HEADER.size:
                fail("STAGE_TRACE_HEADER_FAIL path={}".format(path))
            header = HEADER.unpack(payload)
            magic, kind, sequence, outer_step, grid_level, n_eq, n_arrays, n_items = header
            if magic != MAGIC or kind not in ARRAY_NAMES:
                fail("STAGE_TRACE_HEADER_FAIL path={}".format(path))
            if n_arrays != len(ARRAY_NAMES[kind]) or n_eq == 0 or n_items == 0:
                fail("STAGE_TRACE_SHAPE_FAIL path={}".format(path))
            item_count = n_eq * n_items
            arrays = []
            for _ in range(n_arrays):
                raw = read_exact(stream, item_count * 8, path)
                arrays.append(struct.unpack("<{}d".format(item_count), raw))
            count_records += 1
            yield {
                "key": (kind, sequence, outer_step, grid_level, n_eq, n_items),
                "kind": kind,
                "n_items": n_items,
                "arrays": arrays,
            }
    if count_records == 0:
        fail("STAGE_TRACE_EMPTY path={}".format(path))


def validate_sequence(name, kinds, tolerance, sequence):
    if not all(kind in kinds for kind in (1, 2, 4)):
        fail("STAGE_SEMANTICS_RECORD_FAIL candidate={} sequence={}".format(
            name, sequence))

    face = kinds[1]
    residual = kinds[2]
    metadata = kinds[4]
    n_faces = face["n_items"]
    n_cells = residual["n_items"]
    n_eq = face["key"][4]
    if metadata["n_items"] != n_faces or metadata["key"][4] != 1:
        fail("STAGE_SEMANTICS_SHAPE_FAIL candidate={} sequence={}".format(
            name, sequence))

    left, right, boundary, normal_x, normal_y, normal_z, mesh_vn, area = (
        metadata["arrays"])
    for metadata_name, values in zip(ARRAY_NAMES[4], metadata["arrays"]):
        if not all(math.isfinite(value) for value in values):
            fail("STAGE_SEMANTICS_FINITE_FAIL candidate={} array={}".format(
                name, metadata_name))
    if min(area) <= 0.0:
        fail("STAGE_SEMANTICS_AREA_FAIL candidate={}".format(name))

    left_cells = []
    right_cells = []
    boundary_mask = []
    for face_index in range(n_faces):
        left_cell = int(left[face_index])
        right_cell = int(right[face_index])
        is_boundary = int(boundary[face_index])
        if (left_cell != left[face_index] or right_cell != right[face_index]
                or is_boundary != boundary[face_index]
                or is_boundary not in (0, 1)
                or left_cell < 0 or left_cell >= n_cells
                or (not is_boundary
                    and (right_cell < 0 or right_cell >= n_cells))):
            fail("STAGE_SEMANTICS_CONNECTIVITY_FAIL candidate={} face={}".format(
                name, face_index))
        left_cells.append(left_cell)
        right_cells.append(right_cell)
        boundary_mask.append(is_boundary)

    expected = [0.0] * (n_eq * n_cells)
    invflux = face["arrays"][2]
    for equation in range(n_eq):
        for face_index in range(n_faces):
            value = invflux[equation * n_faces + face_index]
            expected[equation * n_cells + left_cells[face_index]] -= value
            if not boundary_mask[face_index]:
                expected[equation * n_cells + right_cells[face_index]] += value

    actual = residual["arrays"][0]
    residual_error = max(abs(a - b) for a, b in zip(expected, actual))
    conservation_error = 0.0
    for equation in range(n_eq):
        residual_sum = sum(actual[equation * n_cells:(equation + 1) * n_cells])
        boundary_sum = -sum(
            invflux[equation * n_faces + face_index]
            for face_index in range(n_faces)
            if boundary_mask[face_index])
        conservation_error = max(
            conservation_error, abs(residual_sum - boundary_sum))
    if residual_error > tolerance or conservation_error > tolerance:
        fail("STAGE_SEMANTICS_FAIL candidate={} sequence={} residual={:.17g} "
             "conservation={:.17g}".format(
                 name, sequence, residual_error, conservation_error))
    print("STAGE_SEMANTICS_PASS candidate={} sequence={} residual={:.17g} "
          "conservation={:.17g} min_area={:.17g}".format(
              name, sequence, residual_error, conservation_error, min(area)))


def validate_semantics_stream(name, path, tolerance):
    current_sequence = None
    kinds = {}
    record_count = 0
    for record in iter_trace(path):
        sequence = record["key"][1]
        if current_sequence is None:
            current_sequence = sequence
        elif sequence != current_sequence:
            if sequence < current_sequence:
                fail("STAGE_SEMANTICS_ORDER_FAIL candidate={} sequence={}".format(
                    name, sequence))
            validate_sequence(name, kinds, tolerance, current_sequence)
            current_sequence = sequence
            kinds = {}
        if record["kind"] in kinds:
            fail("STAGE_SEMANTICS_DUPLICATE_FAIL candidate={} sequence={} kind={}".format(
                name, sequence, KIND_NAMES[record["kind"]]))
        kinds[record["kind"]] = record
        record_count += 1
    if current_sequence is None:
        fail("STAGE_TRACE_EMPTY path={}".format(path))
    validate_sequence(name, kinds, tolerance, current_sequence)
    return record_count


def compare_record(reference_name, reference, candidate_name, candidate,
                   abs_tol, rel_tol):
    if reference["key"] != candidate["key"]:
        fail("STAGE_TRACE_RECORD_FAIL reference={} candidate={} key_ref={} key_candidate={}".format(
            reference_name, candidate_name, reference["key"], candidate["key"]))

    overall_absolute = 0.0
    overall_scaled = 0.0
    kind = reference["kind"]
    n_items = reference["n_items"]
    for array_name, left, right in zip(
            ARRAY_NAMES[kind], reference["arrays"], candidate["arrays"]):
        if not all(math.isfinite(value) for value in left + right):
            fail("STAGE_TRACE_FINITE_FAIL candidate={} array={}".format(
                candidate_name, array_name))
        absolute = max(abs(a - b) for a, b in zip(left, right))
        scaled = max(
            abs(a - b) / max(abs(a), abs(b), 1.0)
            for a, b in zip(left, right))
        if any(
                abs(a - b) > abs_tol + rel_tol * max(abs(a), abs(b), 1.0)
                for a, b in zip(left, right)):
            fail("STAGE_TRACE_TOLERANCE_FAIL reference={} candidate={} "
                 "kind={} array={} max_absolute={:.17g} max_scaled={:.17g}".format(
                     reference_name, candidate_name, KIND_NAMES[kind], array_name,
                     absolute, scaled))
        overall_absolute = max(overall_absolute, absolute)
        overall_scaled = max(overall_scaled, scaled)
        print("STAGE_TRACE candidate={} kind={} sequence={} array={} n={} "
              "max_absolute={:.17g} max_scaled={:.17g}".format(
                  candidate_name, KIND_NAMES[kind], reference["key"][1],
                  array_name, len(left), absolute, scaled))

        if array_name in ("qf1", "qf2", "state"):
            density = right[:n_items]
            pressure = right[4 * n_items:5 * n_items]
            if min(density) <= 0.0 or min(pressure) <= 0.0:
                fail("STAGE_TRACE_PHYSICAL_FAIL candidate={} array={}".format(
                    candidate_name, array_name))
            print("STAGE_PHYSICAL candidate={} array={} finite=true "
                  "min_density={:.17g} min_pressure={:.17g}".format(
                      candidate_name, array_name, min(density), min(pressure)))
    return overall_absolute, overall_scaled


def compare_streams(reference_name, reference_path, candidates, abs_tol, rel_tol):
    reference = iter_trace(reference_path)
    candidate_iters = [(name, iter_trace(path)) for name, path in candidates]
    count = 0
    aggregates = {name: [0.0, 0.0] for name, _ in candidate_iters}
    while True:
        reference_record = next(reference, None)
        candidate_records = [
            (name, next(iterator, None))
            for name, iterator in candidate_iters
        ]
        if reference_record is None:
            if any(record is not None for _, record in candidate_records):
                fail("STAGE_TRACE_LENGTH_FAIL reference={}".format(reference_name))
            break
        if any(record is None for _, record in candidate_records):
            fail("STAGE_TRACE_LENGTH_FAIL reference={}".format(reference_name))
        for name, candidate_record in candidate_records:
            absolute, scaled = compare_record(
                reference_name, reference_record, name, candidate_record,
                abs_tol, rel_tol)
            aggregates[name][0] = max(aggregates[name][0], absolute)
            aggregates[name][1] = max(aggregates[name][1], scaled)
        count += 1

    for name, _ in candidate_iters:
        print("STAGE_TRACE_PAIR_PASS reference={} candidate={} records={} "
              "max_absolute={:.17g} max_scaled={:.17g}".format(
                  reference_name, name, count,
                  aggregates[name][0], aggregates[name][1]))
    return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("legacy_trace")
    parser.add_argument("cpu_batch_trace")
    parser.add_argument("hip_batch_trace")
    parser.add_argument("--absolute-tolerance", type=float, default=1.0e-12)
    parser.add_argument("--relative-tolerance", type=float, default=1.0e-12)
    args = parser.parse_args()

    semantic_tolerance = max(args.absolute_tolerance, args.relative_tolerance)
    validate_semantics_stream(
        "legacy", args.legacy_trace, semantic_tolerance)
    validate_semantics_stream(
        "cpu_batch", args.cpu_batch_trace, semantic_tolerance)
    validate_semantics_stream(
        "hip_batch", args.hip_batch_trace, semantic_tolerance)
    records = compare_streams(
        "legacy", args.legacy_trace,
        (("cpu_batch", args.cpu_batch_trace),
         ("hip_batch", args.hip_batch_trace)),
        args.absolute_tolerance, args.relative_tolerance)
    print("STAGE_TRACE_PASS records={}".format(records))


if __name__ == "__main__":
    main()
