import sys
import shutil
import importlib


def import_install(module_name):
    """
    Attempt to import a Python module by name; if missing, automatically install it and re-import.

    Parameters
    ----------
    module_name : str
        Name of the module to import and, if necessary, install via pip.

    Behavior
    --------
    - Tries to import the specified module using `importlib`.
    - If the module is not found (`ImportError`), attempts to install it using pip:
      - If running inside a Jupyter environment, uses the `%pip install` magic command.
      - Otherwise, runs `pip install` via a subprocess.
    - After installation, tries to import the module again.
    - Prints status messages indicating success or failure.

    Notes
    -----
    - This function requires internet access and appropriate permissions to install packages.
    - Designed to simplify dependencies management, especially in interactive environments.
    """

    try:
        importlib.import_module(module_name)
        print("success importing ", module_name)
    except ImportError:
        try:
            print(f"Module '{module_name}' not found. Attempting to install...")
            # Use %pip only in Jupyter; fallback otherwise
            import IPython

            ipy = IPython.get_ipython()
            if ipy is not None:
                ipy.run_line_magic("pip", f"install {module_name}")
            else:
                # fallback if not in Jupyter
                import subprocess

                subprocess.check_call(
                    [sys.executable, "-m", "pip", "install", module_name]
                )
            # Try importing again after installation
            importlib.import_module(module_name)
        except Exception as e:
            print(f"<U+274C> Failed to install or import '{module_name}': {e}")


def validate_required_args(required_args, context):
    for argname in required_args:
        if not getattr(context, argname):
            raise ValueError(f"--{argname} must be provided with --{context._action}")


def check_program_exists(exe):
    if shutil.which(exe) is None:
        print(f"Error: '{exe}' not found in PATH.", file=sys.stderr)
        sys.exit(1)


def fit_line(pt1, pt2):
    "return a, b for y=ax+b line through points pt1 = (x1,y1) and pt2 = (x2,y2)"
    dx = pt1[0] - pt2[0]
    dy = pt1[1] - pt2[1]
    # special case: vertical line
    if dx == 0:
        print("Vertical cut at ", pt1[0])
        return (None, None)
    # sloping line
    slope = dy / dx
    crosspoint = pt2[1] - slope * pt2[0]
    return (slope, crosspoint)
