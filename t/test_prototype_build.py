#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class PrototypeBuildTests(unittest.TestCase):
    def test_standalone_build_refreshes_libraries_once_and_only_relinks_changes(self):
        with tempfile.TemporaryDirectory(prefix="starfix-build-") as directory:
            root = Path(directory)
            prototype = root / "prototype"
            prototype.mkdir()
            shutil.copyfile(ROOT / "prototype/Makefile", prototype / "Makefile")
            (root / "Makefile").write_text(
                ".PHONY: record\n"
                "record:\n\t@echo checked >> checks.log\n"
                "libstarfix.a: core.c | record\n\tcp $< $@\n"
                "libstarfix_navigation.a: navigation.c | record\n\tcp $< $@\n"
            )
            for name in ("core.c", "navigation.c"):
                (root / name).write_text(name)
            for name in ("main", "camera", "bno055", "sim", "camera_mock", "bno055_mock",
                         "trail_preprocessor"):
                (prototype / f"{name}.c").write_text(name)
            compiler = root / "compiler.py"
            compiler.write_text(
                "import pathlib, sys\n"
                "output = pathlib.Path(sys.argv[sys.argv.index('-o') + 1])\n"
                "output.write_text('built')\n"
                "with open('compiler.log', 'a') as log:\n"
                "    log.write(output.name + '\\n')\n"
            )

            def build():
                # The nested make is independent of any outer make check jobserver.
                env = os.environ.copy()
                env.pop("MAKEFLAGS", None)
                env.pop("MFLAGS", None)
                result = subprocess.run(
                    ["make", "-j8", "-C", str(prototype),
                     f"CC={sys.executable} {compiler}"],
                    capture_output=True, text=True, env=env, timeout=30,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            def log(path):
                return path.read_text().splitlines() if path.exists() else []

            checks = root / "checks.log"
            compilations = prototype / "compiler.log"
            build()
            self.assertEqual(log(checks), ["checked"])
            self.assertIn("starfix_prototype", log(compilations))
            self.assertIn("starfix_prototype_sim", log(compilations))

            # Fixed old timestamps make the update checks independent of clock resolution.
            baseline = time.time() - 3600
            for path in root.rglob("*"):
                if path.is_file():
                    os.utime(path, (baseline, baseline))
            before = log(compilations)
            build()
            self.assertEqual(log(checks), ["checked"] * 2)
            self.assertEqual(log(compilations), before, "an unchanged build must not relink")

            for index, (source, archive) in enumerate(
                (("core.c", "libstarfix.a"),
                 ("navigation.c", "libstarfix_navigation.a")), start=3,
            ):
                with self.subTest(source=source):
                    # Leave objects/binaries older than the changed root source.
                    for path in root.rglob("*"):
                        if path.is_file():
                            os.utime(path, (baseline, baseline))
                    changed = root / source
                    changed.write_text("updated " + source)
                    os.utime(changed, (baseline + 1, baseline + 1))
                    before = log(compilations)
                    build()
                    self.assertEqual((root / archive).read_text(), changed.read_text())
                    self.assertEqual(log(checks), ["checked"] * index)
                    self.assertCountEqual(log(compilations)[len(before):],
                                          ["starfix_prototype", "starfix_prototype_sim"])

            before = log(compilations)
            build()
            self.assertEqual(log(checks), ["checked"] * 5)
            self.assertEqual(log(compilations), before)


if __name__ == "__main__":
    unittest.main()
