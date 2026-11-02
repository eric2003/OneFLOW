"""Serialize validated layout topology into OneFLOW's legacy layout syntax."""


def serialize_topology(topology, format_coordinate, format_distribution):
    """Render shared topology records while delegating distribution syntax."""
    output = []

    for point in topology.points:
        coordinates = ",".join(
            format_coordinate(value) for value in (point.x, point.y, point.z)
        )
        output.append(f"Point({point.id})= {{{coordinates}}}")

    output.append("")
    for line in topology.lines:
        output.append(f"Line({line.id})={{{line.start_point},{line.end_point}}}")

    output.append("")
    for line in topology.lines:
        output.append(f"Dim({line.id})={{{line.point_count}}}")

    output.append("")
    for line in topology.lines:
        output.append(format_distribution(line))

    output.append("")
    for line in topology.lines:
        output.append(f"Boundary({line.id})={{{line.boundary_type}}}")

    for face in topology.faces:
        output.append("")
        for position, line_id in face.line_by_side:
            output.append(f"Add L2F {face.id} {position} {line_id}")

    output.append("")
    for block in topology.blocks:
        for position, face_id in block.face_by_position:
            output.append(f"Add F2B {block.id} {position} {face_id}")

    output.append("")
    return "\n".join(output)


__all__ = ["serialize_topology"]
