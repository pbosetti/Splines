#!/usr/bin/env python3
"""
Smoke test for the Splines Python wrapper.

Usage, from the Splines repository root after building/installing the wrapper:

    conda activate splines-py
    python test_splines.py

If the module is not on PYTHONPATH, either run:

    export PYTHONPATH="$PWD/lib/lib:$PYTHONPATH"
    python test_splines.py

or pass the directory containing Splines*.so:

    python test_splines.py --module-dir ./lib/lib

The script is deliberately defensive: pybind11 wrappers for this project may
expose slightly different method names depending on the branch/version.
"""

from __future__ import annotations

import argparse
import importlib
import inspect
import math
import pathlib
import sys
import traceback
from typing import Any, Iterable


def add_module_dir(path: str | None) -> None:
    if path:
        p = pathlib.Path(path).expanduser().resolve()
        sys.path.insert(0, str(p))
        print(f"[info] added to sys.path: {p}")


def import_splines() -> Any:
    print("[test] importing Splines")
    mod = importlib.import_module("Splines")
    print(f"[ok] imported Splines from: {getattr(mod, '__file__', '<built-in>')}")
    return mod


def public_names(obj: Any) -> list[str]:
    return sorted(name for name in dir(obj) if not name.startswith("_"))


def show_public_api(mod: Any) -> None:
    names = public_names(mod)
    print("\n[info] public symbols in Splines:")
    for name in names:
        attr = getattr(mod, name)
        kind = "class" if inspect.isclass(attr) else "callable" if callable(attr) else type(attr).__name__
        print(f"  - {name} ({kind})")


def try_call(description: str, fn, *args, **kwargs) -> tuple[bool, Any]:
    print(f"\n[test] {description}")
    try:
        out = fn(*args, **kwargs)
    except Exception as exc:  # noqa: BLE001 - smoke-test script
        print(f"[fail] {description}: {type(exc).__name__}: {exc}")
        return False, exc
    print(f"[ok] {description}")
    return True, out


def first_existing_method(obj: Any, candidates: Iterable[str]) -> str | None:
    for name in candidates:
        if hasattr(obj, name) and callable(getattr(obj, name)):
            return name
    return None


def try_build_spline_class(cls: type, x: list[float], y: list[float]) -> tuple[bool, Any]:
    """Try common construction patterns used by C++/pybind spline wrappers."""
    class_name = cls.__name__

    # Pattern 1: constructor with no arguments, then build/setup method.
    ok, obj_or_exc = try_call(f"construct {class_name}()", cls)
    if ok:
        obj = obj_or_exc
        print(f"[info] methods on {class_name}: {public_names(obj)}")
        method_name = first_existing_method(
            obj,
            [
                "build",
                "setup",
                "setData",
                "set_data",
                "buildSpline",
                "build_spline",
                "init",
            ],
        )
        if method_name is not None:
            method = getattr(obj, method_name)
            for args in ((x, y), (len(x), x, y)):
                ok2, _ = try_call(f"{class_name}.{method_name}{args!r}", method, *args)
                if ok2:
                    return True, obj
        else:
            print(f"[warn] no known build/setup method found on {class_name}")

    # Pattern 2: constructor directly from x, y.
    ok, obj_or_exc = try_call(f"construct {class_name}(x, y)", cls, x, y)
    if ok:
        return True, obj_or_exc

    # Pattern 3: constructor from n, x, y.
    ok, obj_or_exc = try_call(f"construct {class_name}(n, x, y)", cls, len(x), x, y)
    if ok:
        return True, obj_or_exc

    return False, None


def try_evaluate_spline(obj: Any, points: list[float]) -> bool:
    eval_name = first_existing_method(
        obj,
        [
            "eval",
            "evaluate",
            "operator",
            "value",
            "Y",
            "y",
            "__call__",
        ],
    )
    if eval_name is None:
        print("[warn] no known evaluation method found")
        return False

    method = getattr(obj, eval_name)
    ok_all = True
    for xx in points:
        ok, value = try_call(f"evaluate at x={xx} using {eval_name}", method, xx)
        if ok:
            print(f"       {eval_name}({xx}) = {value}")
            if isinstance(value, (float, int)) and not math.isfinite(float(value)):
                print("[fail] non-finite value")
                ok_all = False
        else:
            ok_all = False
    return ok_all


def smoke_test_classes(mod: Any) -> bool:
    # Simple monotone test data.
    x = [0.0, 1.0, 2.0, 3.0]
    y = [0.0, 1.0, 4.0, 9.0]
    points = [0.0, 0.5, 1.5, 3.0]

    preferred = [
        "LinearSpline",
        "CubicSpline",
        "AkimaSpline",
        "VanLeerSpline",
        "PchipSpline",
        "QuinticSpline",
        "ConstantSpline",
    ]

    found_any_class = False
    built_any = False

    for name in preferred:
        cls = getattr(mod, name, None)
        if cls is None or not inspect.isclass(cls):
            continue
        found_any_class = True
        print(f"\n========== Testing class {name} ==========")
        ok, obj = try_build_spline_class(cls, x, y)
        if ok:
            built_any = True
            try_evaluate_spline(obj, points)

    if not found_any_class:
        print("[warn] no expected spline classes found by name")
    if not built_any:
        print("[warn] no spline object could be built with the generic patterns")

    return found_any_class


def main() -> int:
    parser = argparse.ArgumentParser(description="Smoke test for the Splines Python wrapper")
    parser.add_argument(
        "--module-dir",
        default=None,
        help="directory containing Splines*.so, e.g. ./lib/lib",
    )
    parser.add_argument(
        "--api-only",
        action="store_true",
        help="only import the module and print the public API",
    )
    args = parser.parse_args()

    add_module_dir(args.module_dir)

    try:
        mod = import_splines()
        show_public_api(mod)
        if not args.api_only:
            smoke_test_classes(mod)
    except Exception:  # noqa: BLE001 - command-line smoke test
        print("\n[error] smoke test failed with an unexpected exception:")
        traceback.print_exc()
        return 1

    print("\n[done] smoke test completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
