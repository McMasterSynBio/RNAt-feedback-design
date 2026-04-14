# RNAt-feedback-design
Feedback design system to optimize RNA thermosensor sequences for desired level of stability and controlled melting.

- NUPACK Package: https://www.nupack.org/download/software 
- NUPACK Install: https://docs.nupack.org/start/
- uv Install: https://docs.astral.sh/uv/guides/projects/#uvlock
- ViennaRNA Install: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/install.html
- nuad Install: https://nuad.readthedocs.io/en/latest/

Below are the instructions on how to set up the uv project:

1. For the installation of `uv`, please refer to [the installation page](https://docs.astral.sh/uv/getting-started/installation/)

2. After navigating to the project repository root, run `uv sync` to set up a virtual environment and install the dependencies.

3. If you're on VS Code, please ensure that your IDE recognizes the python directory inside the virtual environment (under `.venv`) to run the scripts.

4. You can now run the main script. An example instance `uv run python main.py --target-tm 42 --seq-length 25 --num-runs 5`
