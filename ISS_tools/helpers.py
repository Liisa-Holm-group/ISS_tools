import sys
import shutil


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
