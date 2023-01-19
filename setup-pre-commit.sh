# Set up pre-commit hooks, in this case just clang-format checking
#
# Note: this overwrites existing .git/hooks/pre-commit
#
# See .clang-format for configuration file
# Code copied from arg_needle_lib,
# modified from https://github.com/KDAB/kdabtv/tree/master/Qt-Widgets-and-more/clang-format
# Based on this tutorial: https://www.youtube.com/watch?v=Cz36YveDI2E

echo "#!/bin/sh
python3 .git/hooks/pre-commit-clang-format.py" > .git/hooks/pre-commit


echo "import subprocess
try:
    output = str(subprocess.check_output([\"git\", \"clang-format\", \"--diff\"]))
except subprocess.CalledProcessError as e:
    print(e)
    print(\"Error raised, try installing clang-format.\\n\")
    exit(1)
 
if \"clang-format did not modify any files\" not in output and \"no modified files to format\" not in output:
    print(\"Run git clang-format, add the modified files, then commit.\\n\")
    exit(1)
else:
    exit(0)" > .git/hooks/pre-commit-clang-format.py


chmod +x .git/hooks/pre-commit .git/hooks/pre-commit-clang-format.py