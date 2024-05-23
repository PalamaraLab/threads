# Building and debugging in VSCode devcontainer

The devcontainer folder provides a consistent build and debug environment that can be run in VSCode locally or in GitHub Codespaces remotely.

## Opening devcontainer locally

When you open the cloned TDPBWT folder in VSCode, accept the popup to open in a Dev Container or run the editor command 'Dev Containers: Reopen in Dev Container'.

## Initial setup

Threads has python dependencies and the simplest way to get these is to simply build everything with:

```sh
pip install .
```

Along with getting the depenencies, this builds the Threads binary into `/usr/local/lib/python3.10/dist-packages/`, which has the side-effect of setting up CMake to build binaries in `Release` format (see `CMAKE_BUILD_TYPE` in `build/CMakeCache.txt`).

In order to debug the code, In VSCode, run the editor command `CMake: Select Variant` and select either `Debug` or `RelWithDebInfo`.

## Building and running

In VSCode, run the editor command 'CMake: Build' to rebuild into `build/src`. To use this version in python, set the `PYTHONPATH` environment variable to `build/src` so this folder is searched for modules before `dist-packages`. Running without this environment variable will fall back to the original pip-installed version.

The unit tests launch the build C++ directly. The python launch works slightly differently, by running the python interpreter (with `PYTHONPATH` set correctly) so the debug symbols available are when the `.so` file is imported by python; the debugger will stop on any breakpoints you have set.

## Running outside of a devcontainer

This build and debug configuration also works fine outside of a devcontainer provided the correct dependencies are installed (see `Dockerfile` for reference).

## `.vscode` folder

This folder contains everything to immediately run and debug in VSCode.

- `settings.json` shows common options for Debugging that can be uncommented to enable.
- `launch.json` contains two launch configurations:
  - Run/debug unit tests.
  - Run/debug infer example in main README.

## Optionally building arg-needle-lib

If you are running on a platform where wheels do not exist for arg-needle-lib, clone it and build the wheel locally with `pip install .` in the source folder (which may also require `apt install libboost-iostreams-dev libgsl-dev`). Then `pip install` the generated `arg_needle_lib*.whl` so it can be picked up by threads.

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