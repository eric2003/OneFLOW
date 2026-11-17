"""Tests for lossless legacy parameter-file parsing and updates."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.parameter_document import ParameterDocument, parameter_specs_for_file


class ParameterDocumentTest(unittest.TestCase):
    def test_reads_task_selector_and_specs(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "simu.txt"
            path.write_bytes(b'string simutask = "Solve";\r\n')
            document = ParameterDocument.load(path)
            self.assertEqual(document.get("simutask").value, "Solve")
            self.assertEqual(parameter_specs_for_file(path)[0].key, "simutask")

    def test_gbk_update_preserves_chinese_comments_and_crlf(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "cfd.txt"
            comments = (
                "\u6700\u5927\u8fed\u4ee3\u6b65\u6570",
                "\u6700\u5927\u65f6\u95f4",
                "\u4fdd\u7559",
            )
            original = (
                "int maxSteps = 50; //" + comments[0] + "\r\n"
                "Real maxTime = 0.2; //" + comments[1] + "\r\n"
                "int ischeme = 1; //Roe\r\n"
                "int unknownOption = 7; //" + comments[2] + "\r\n"
            )
            path.write_bytes(original.encode("gbk"))
            document = ParameterDocument.load(path)
            self.assertEqual(document.encoding, "gbk")
            self.assertTrue(document.save({"maxSteps": 120, "ischeme": 12}))
            updated = path.read_bytes().decode("gbk")
            self.assertIn("int maxSteps = 120; //" + comments[0] + "\r\n", updated)
            self.assertIn("int ischeme = 12; //Roe\r\n", updated)
            self.assertIn("int unknownOption = 7; //" + comments[2] + "\r\n", updated)
            self.assertEqual(updated.count("\r\n"), 4)

    def test_utf8_bom_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "simu.txt"
            path.write_bytes(b'\xef\xbb\xbfstring simutask = "Grid";\n')
            document = ParameterDocument.load(path)
            self.assertEqual(document.encoding, "utf-8-sig")
            self.assertTrue(document.save({"simutask": "Solve"}))
            self.assertTrue(path.read_bytes().startswith(b"\xef\xbb\xbf"))
            self.assertEqual(path.read_bytes().decode("utf-8-sig"), 'string simutask = "Solve";\n')

    def test_unknown_text_comments_and_crlf_are_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "cfd.txt"
            original = b"// header\r\n\r\nint maxSteps = 50; // count\r\ncustom syntax\r\n"
            path.write_bytes(original)
            document = ParameterDocument.load(path)
            self.assertTrue(document.save({"maxSteps": 75}))
            with path.open("rb") as stream:
                updated = stream.read()
            self.assertEqual(updated, b"// header\r\n\r\nint maxSteps = 75; // count\r\ncustom syntax\r\n")

    def test_rejects_concurrent_external_edit(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "simu.txt"
            path.write_text('string simutask = "Grid";\n', encoding="utf-8")
            document = ParameterDocument.load(path)
            path.write_text('string simutask = "Solve";\n', encoding="utf-8")
            with self.assertRaises(RuntimeError):
                document.save({"simutask": "Grid"})

    def test_noop_does_not_reformat(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "cfd.txt"
            original = b"int maxSteps   =  50; // spacing\n"
            path.write_bytes(original)
            document = ParameterDocument.load(path)
            self.assertFalse(document.save({"maxSteps": 50}))
            self.assertEqual(path.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
