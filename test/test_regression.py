import subprocess
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent


def run_threads_infer():

    command = [
        "threads", "infer",
        "--pgen", str(BASE_DIR / "example/example_data.pgen"),
        "--map_gz", str(BASE_DIR / "example/example_data.map"),
        "--demography", str(BASE_DIR / "example/Ne10000.demo"),
        "--out", str(BASE_DIR / "example/example_data.threads")
    ]

    return subprocess.run(command, cwd=BASE_DIR, capture_output=True, text=True)


def run_threads_convert():

    command = [
        "threads", "convert",
        "--threads", str(BASE_DIR / "example/example_data.threads"),
        "--argn", str(BASE_DIR / "example/example_data.argn")
    ]

    return subprocess.run(command, cwd=BASE_DIR, capture_output=True, text=True)


def test_threads_on_example_data():

    result_infer = run_threads_infer()
    assert result_infer.returncode == 0, f"Script did not run successfully: {result_infer.stderr}"
    assert (BASE_DIR / 'example' / 'example_data.threads').is_file()

    result_convert = run_threads_convert()
    assert result_convert.returncode == 0, f"Script did not run successfully: {result_convert.stderr}"
    assert (BASE_DIR / 'example' / 'example_data.argn').is_file()
