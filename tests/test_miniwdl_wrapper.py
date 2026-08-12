from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
WRAPPER = REPOSITORY_ROOT / "deploy" / "miniwdl" / "miniwdl"


FAKE_MINIWDL = r'''#!/usr/bin/env python3
import json
import os
import sys
from pathlib import Path

state = Path(os.environ["FAKE_MINIWDL_STATE"])
attempt = int(state.read_text()) + 1 if state.exists() else 1
state.write_text(str(attempt))
mode = os.environ["FAKE_MINIWDL_MODE"]
output = Path(sys.argv[sys.argv.index("-o") + 1])

race = """Traceback (most recent call last):
  File "/opt/lib/python/site-packages/WDL/runtime/cache.py", line 1, in get
    self._check_files_coherence(cache_file)
  File "/opt/lib/python/site-packages/WDL/runtime/cache.py", line 2, in _check_files_coherence
    self._effective_mtime(cache_file)
  File "/opt/lib/python/site-packages/WDL/runtime/cache.py", line 3, in _effective_mtime
    os.stat(cache_file)
FileNotFoundError: [Errno 2] No such file or directory: '/hpc/temp/user/.miniwdl/miniwdl_call_cache/task/key.json'
"""

if mode == "race_then_success" and attempt == 1:
    print(race, end="")
    sys.exit(2)
if mode in ("race_then_success", "success"):
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"outputs": {"result": "ok"}}))
    print("workflow complete")
    sys.exit(0)
if mode == "always_race":
    print(race, end="")
    sys.exit(2)
if mode == "missing_output":
    print("workflow claimed success")
    sys.exit(0)
print("task failed in container")
sys.exit(42)
'''


class MiniwdlWrapperTest(unittest.TestCase):
    def _run(self, root: Path, mode: str) -> subprocess.CompletedProcess[str]:
        install_root = root / "miniwdl"
        bin_directory = install_root / "bin"
        env_directory = install_root / "env" / "bin"
        config_directory = install_root / "config"
        bin_directory.mkdir(parents=True)
        env_directory.mkdir(parents=True)
        config_directory.mkdir(parents=True)
        wrapper = bin_directory / "miniwdl"
        shutil.copy2(WRAPPER, wrapper)
        wrapper.chmod(0o755)
        fake = env_directory / "miniwdl"
        fake.write_text(textwrap.dedent(FAKE_MINIWDL), encoding="utf-8")
        fake.chmod(0o755)
        os.symlink(shutil.which("python3"), env_directory / "python")
        (config_directory / "miniwdl.cfg").write_text("", encoding="utf-8")

        wdl = root / "workflow.wdl"
        inputs = root / "input.json"
        output = root / "outputs.json"
        run_dir = root / "runs"
        wdl.write_text("version 1.0\n", encoding="utf-8")
        inputs.write_text("{}\n", encoding="utf-8")
        environment = os.environ.copy()
        environment.update(
            {
                "FAKE_MINIWDL_MODE": mode,
                "FAKE_MINIWDL_STATE": str(root / "attempt-count"),
                "LONGPIPE_MINIWDL_BACKOFF_SECONDS": "0",
                "MINIWDL_CACHE_ROOT": str(root / "cache"),
            }
        )
        return subprocess.run(
            [
                str(wrapper),
                "run",
                str(wdl),
                "-i",
                str(inputs),
                "--dir",
                str(run_dir),
                "-o",
                str(output),
                "--verbose",
            ],
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )

    def test_cache_race_retries_then_succeeds(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            completed = self._run(root, "race_then_success")
            output = json.loads((root / "outputs.json").read_text())

            self.assertEqual(completed.returncode, 0, completed.stdout)
            self.assertEqual((root / "attempt-count").read_text(), "2")
            self.assertEqual(output["outputs"]["result"], "ok")
            self.assertTrue(
                (root / "runs" / "longpipe-miniwdl-attempt-1.log").is_file()
            )
            self.assertTrue(
                (root / "runs" / "longpipe-miniwdl-attempt-2.log").is_file()
            )
            self.assertIn(
                "Detected transient MiniWDL call-cache race",
                completed.stdout,
            )

    def test_nonmatching_failure_is_not_retried(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            completed = self._run(root, "ordinary_failure")

            self.assertEqual(completed.returncode, 42)
            self.assertEqual((root / "attempt-count").read_text(), "1")
            self.assertFalse(
                (root / "runs" / "longpipe-miniwdl-attempt-2.log").exists()
            )
            self.assertIn("does not match", completed.stdout)

    def test_repeated_cache_race_exhausts_retry(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            completed = self._run(root, "always_race")

            self.assertEqual(completed.returncode, 2)
            self.assertEqual((root / "attempt-count").read_text(), "2")
            self.assertTrue(
                (root / "runs" / "longpipe-miniwdl-attempt-2.log").is_file()
            )
            self.assertIn("persisted through all attempts", completed.stdout)

    def test_success_requires_valid_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            completed = self._run(root, "missing_output")

            self.assertEqual(completed.returncode, 1)
            self.assertEqual((root / "attempt-count").read_text(), "1")
            self.assertIn("output is missing or empty", completed.stdout)


if __name__ == "__main__":
    unittest.main()
