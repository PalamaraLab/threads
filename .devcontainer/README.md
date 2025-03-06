# Building and debugging in VSCode devcontainer

The devcontainer folder provides a consistent build and debug environment that can be run in VSCode locally or in GitHub Codespaces remotely.

## Local devcontainer

After cloning the repo, open in VSCode and accept the popup to open in a Dev Container. You can also run this using the editor command 'Dev Containers: Reopen in Dev Container'.

## Remote devcontainer

Alternatively, you can run and debug remotely in GitHub by starting a new Codespace. In the top left menu, select Codespaces and click on New Codespace. Then select this repo name and the branch you are working on, then Create Codespace. This will spin up a browser-based VSCode instance running on a remote machine.

## Python setup

The simplest approach is to install from the local folder with:

```sh
pip install -e .[dev]
```

This will fetch dependencies, compile this project and install them. The `[dev]` option ensures that `pytest` is available for running unit tests.

The `-e` editable mode flag ensures that changes to local files in the project are picked up by unit tests without having to run `pip install` again.

## C++ setup

Follow the Python setup above to install all dependencies, then open in VSCode and open in a Dev Container as above.

Run the VSCode editor command 'CMake: Build' which will configure CMake to build with Debug information and build the binary into `build/src`. If prompted to select the right kit to build with then select the default.

Because we have run `pip install .` above, by default any tests will try to run that install rather than our locally-built binaries. There are two different ways to make it work with our code changes:

1. Set the `PYTHONPATH` environment variable to `src:build/src`. This will ensure that scripts that `import threads_arg` will look for Python files in `src/threads_arg` and and the corresponding `.so` module binary built by CMake in `build/src`. These folders are searched before `dist-packages`. Running without this environment variable will fall back to the original pip-installed version.
2. To avoid confusion, when working in C++ you may wish to `pip uninstall threads-arg` entirely and always depend on `PYTHONPATH` being set.

If you look in `launch.json` you will see the tests set `PYTHONPATH` as above to ensure the correct files are loaded.

## Running outside of a devcontainer

This build and debug configuration also works fine outside of a devcontainer provided the correct dependencies are installed (see `Dockerfile` for reference).

## Files in `.vscode` folder

These files set up VSCode run and debug options, including example scenarios for debugging C++ unit tests, Python tests, and launching threads from Python whilst debugging the C++ library.

- `settings.json`: Shows common options for Debugging that can be uncommented to enable.
- `launch.json`: Example run and debug scenarios.

## Files in `.devcontainer` folder

These define requirements to spin up a devcontainer

- `devcontainer.json`: Name, VSCode extensions required and system flags for debugging.
- `Dockerfile`: Ubuntu image with all required system and python libraries.

## Additional debug information

When building in GCC, the values of C++ `std` types (for example `ostringstream`) are not displayed correctly. This can be improved by installing the `libstd++` debug libs and a custom gdb printer. In bash:

```sh
apt-get install -y libstdc++6-11-dbg
```

Download [`printer.py`](https://github.com/gcc-mirror/gcc/blob/master/libstdc%2B%2B-v3/python/libstdcxx/v6/printers.py) and add it to your `.gdbinit` file, e.g.

```
python
import sys
sys.path.insert(0, '<parent_dir_for_printer.py>')
from printers import register_libstdcxx_printers
register_libstdcxx_printers (None)
end
```

Additionally, invoke `.gdbinit` before other `launch.json` settings:

```json
    "setupCommands": [
        {
            "text": "source <parent_dir_for_gdbinit>/.gdbinit",
            "ignoreFailures": false
        },
        ...
```

Then `LD_PRELOAD` the debug libs in the `launch.json` environment:

```json
    "environment": [
        {
            "name": "LD_PRELOAD",
            "value": "/usr/lib/x86_64-linux-gnu/debug/libstdc++.so"
        }
    ]
```
