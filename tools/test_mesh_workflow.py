import unittest

from mesh_workflow import validate_summary


class MeshWorkflowValidationTest(unittest.TestCase):
    def setUp(self):
        self.expected = {
            "zones": 1,
            "dimension": 2,
            "nodes": 9,
            "cells": 4,
            "faces": 12,
            "boundaryFaces": 8,
            "internalFaces": 4,
            "boundingBox": {
                "min": [0.0, 0.0, 0.0],
                "max": [1.0, 1.0, 0.0],
            },
        }

    def test_validates_runtime_summary_shape(self):
        summary = {
            "dimension": 2,
            "zones": [{"name": "zone0"}],
            "counts": {
                "nodes": 9,
                "cells": 4,
                "faces": 12,
                "boundaryFaces": 8,
            },
            "boundingBox": {
                "min": [0.0, 0.0, 0.0],
                "max": [1.0, 1.0, 0.0],
            },
        }

        validate_summary(summary, self.expected, "runtime")

    def test_validates_export_summary_shape(self):
        summary = {
            "dimension": 2,
            "counts": {
                "zones": 1,
                "nodes": 9,
                "cells": 4,
                "faces": 12,
                "boundaryFaces": 8,
            },
            "boundingBox": {
                "min": [0.0, 0.0, 0.0],
                "max": [1.0, 1.0, 0.0],
            },
        }

        validate_summary(summary, self.expected, "export")

    def test_rejects_wrong_face_count(self):
        summary = {
            "dimension": 2,
            "counts": {
                "zones": 1,
                "nodes": 9,
                "cells": 4,
                "faces": 11,
                "boundaryFaces": 8,
            },
            "boundingBox": {
                "min": [0.0, 0.0, 0.0],
                "max": [1.0, 1.0, 0.0],
            },
        }

        with self.assertRaisesRegex(ValueError, "export faces"):
            validate_summary(summary, self.expected, "export")

    def test_rejects_wrong_bounding_box(self):
        summary = {
            "dimension": 2,
            "counts": {
                "zones": 1,
                "nodes": 9,
                "cells": 4,
                "faces": 12,
                "boundaryFaces": 8,
            },
            "boundingBox": {
                "min": [0.0, 0.0, 0.0],
                "max": [1.1, 1.0, 0.0],
            },
        }

        with self.assertRaisesRegex(ValueError, "boundingBox.max"):
            validate_summary(summary, self.expected, "export")


if __name__ == "__main__":
    unittest.main()
