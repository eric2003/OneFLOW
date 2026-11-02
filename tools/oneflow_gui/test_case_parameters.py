import tempfile
import unittest
from pathlib import Path
from tools.oneflow_gui.case_parameters import ParameterDocument, discover_parameter_documents

class CaseParameterTest(unittest.TestCase):
    def test_gb18030_round_trip_preserves_comments_and_crlf(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cfd.txt"
            path.write_bytes("int maxSteps = 50; // 总模拟步数\r\n".encode("gb18030"))
            doc = ParameterDocument.read(path)
            self.assertEqual(doc.encoding, "gb18030")
            self.assertEqual(doc.replace_value("maxSteps", 120),
                "int maxSteps = 120; // 总模拟步数\r\n".encode("gb18030"))

    def test_duplicate_declarations_remain_addressable(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cfd.txt"
            path.write_text("int n = 1;\nint n = 2;\n", encoding="utf-8")
            doc = ParameterDocument.read(path)
            self.assertEqual(len(doc.values("n")), 2)
            self.assertIn(b"int n = 2;", doc.replace_value("n", 3))

    def test_manifest_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "control.txt").write_text('"simu.txt";\n"cfd.txt";\n', encoding="ascii")
            (root / "simu.txt").write_text('string simutask = "Solve";\n', encoding="ascii")
            (root / "cfd.txt").write_text("int maxSteps = 50;\n", encoding="ascii")
            self.assertEqual([d.path.name for d in discover_parameter_documents(root / "control.txt")],
                             ["simu.txt", "cfd.txt"])
