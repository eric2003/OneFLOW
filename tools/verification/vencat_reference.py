"""Generate the historical Vencat analytic reference profile."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable


def vencat(x: float) -> float:
    """Return the original Vencat profile with its fixed coefficient."""
    return (1.0 + 2.0 * x) / (1.0 + x + 2.0 * x * x)


def vencat_c(x: float, c: float) -> float:
    """Return the coefficient-parameterized Vencat profile."""
    return (1.0 + c * x) / (1.0 + x + c * x * x)


def vencat_eric(x: float) -> float:
    """Return the historical Eric variant of the Vencat profile."""
    return (1.0 + 3.0 * x) / (1.0 + x + 2.0 * x * x)


def write_vencat_reference(
    output_path: str | Path,
    *,
    point_count: int = 101,
    c_values: Iterable[float] = range(2, 12),
) -> Path:
    """Write one zone per coefficient to an explicitly selected file."""
    if not isinstance(point_count, int) or isinstance(point_count, bool):
        raise TypeError("point_count must be an integer")
    if point_count < 2:
        raise ValueError("point_count must be at least 2")

    coefficients = tuple(float(value) for value in c_values)
    if not coefficients:
        raise ValueError("c_values must contain at least one coefficient")
    if any(not math.isfinite(value) or value <= 0.0 for value in coefficients):
        raise ValueError("all c_values must be finite and positive")

    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)

    with target.open("w", encoding="ascii", newline="\n") as stream:
        stream.write('title="THE FLOW FIELD OF ONEFLOW"\n')
        stream.write("variables=\n")
        stream.write('"x"\n')
        stream.write('"vencat"\n')
        stream.write('"gvencat"\n')

        for coefficient in coefficients:
            stream.write(f" zone  i = {point_count}\n")
            for index in range(point_count):
                x = 2.0 * index / (point_count - 1)
                value = vencat_c(x, coefficient)
                stream.write(f"{x:.17g} {value:.17g} {x * value:.17g}\n")

    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate a multi-zone Vencat analytic reference profile."
    )
    parser.add_argument("output", type=Path, help="explicit output file path")
    parser.add_argument("--points", type=int, default=101)
    parser.add_argument("--c-start", type=float, default=2.0)
    parser.add_argument("--c-count", type=int, default=10)
    args = parser.parse_args(argv)

    if args.c_count < 1:
        parser.error("--c-count must be positive")

    coefficients = (
        args.c_start + offset for offset in range(args.c_count)
    )
    output = write_vencat_reference(
        args.output,
        point_count=args.points,
        c_values=coefficients,
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
