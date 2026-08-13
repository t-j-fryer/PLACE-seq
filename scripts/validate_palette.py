#!/usr/bin/env python3
"""Check a categorical palette for print and colour-vision-deficiency safety.

Figure colour is computable, so it is computed here rather than judged by eye.
Each candidate colour is converted to OKLab, simulated under protanopia and
deuteranopia with the Machado-Oliveira-Fernandes (2009) model at severity 1.0,
and scored on:

1. lightness band       OKLCH L within the mode's band
2. chroma floor         OKLCH C >= 0.10, below which a hue reads as grey
3. CVD separation       pairwise OKLab dE >= 8 (floor 6) under both simulations
4. normal-vision floor  worst pair dE >= 15 unsimulated
5. contrast vs surface  >= 3:1 for marks

``dE`` throughout is Euclidean distance in OKLab scaled by 100.

    python scripts/validate_palette.py "#0072B2,#E69F00" --mode light --pairs all
"""

from __future__ import annotations

import argparse
import itertools
import math
import sys

# Machado, Oliveira & Fernandes (2009), severity 1.0, applied in linear RGB.
_CVD_MATRICES = {
    "protanopia": (
        (0.152286, 1.052583, -0.204868),
        (0.114503, 0.786281, 0.099216),
        (-0.003882, -0.048116, 1.051998),
    ),
    "deuteranopia": (
        (0.367322, 0.860646, -0.227968),
        (0.280085, 0.672501, 0.047413),
        (-0.011820, 0.042940, 0.968881),
    ),
}

LIGHTNESS_BAND = {"light": (0.43, 0.77), "dark": (0.48, 0.67)}
DEFAULT_SURFACE = {"light": "#FFFFFF", "dark": "#101418"}
CHROMA_FLOOR = 0.10
CVD_TARGET, CVD_FLOOR = 8.0, 6.0
NORMAL_VISION_FLOOR = 15.0
CONTRAST_FLOOR = 3.0


def parse_hex(value: str) -> tuple[float, float, float]:
    text = value.strip().lstrip("#")
    if len(text) == 3:
        text = "".join(ch * 2 for ch in text)
    if len(text) != 6:
        raise ValueError(f"not a hex colour: {value!r}")
    return tuple(int(text[i : i + 2], 16) / 255 for i in (0, 2, 4))  # type: ignore[return-value]


def _to_linear(channel: float) -> float:
    return channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4


def _to_srgb(channel: float) -> float:
    value = 12.92 * channel if channel <= 0.0031308 else 1.055 * channel ** (1 / 2.4) - 0.055
    return min(1.0, max(0.0, value))


def linear_rgb(hex_colour: str) -> tuple[float, float, float]:
    return tuple(_to_linear(c) for c in parse_hex(hex_colour))  # type: ignore[return-value]


def oklab(linear: tuple[float, float, float]) -> tuple[float, float, float]:
    r, g, b = linear
    long_ = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b
    medium = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b
    short = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b
    l_, m_, s_ = (math.copysign(abs(v) ** (1 / 3), v) for v in (long_, medium, short))
    return (
        0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_,
        1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_,
        0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_,
    )


def simulate(linear: tuple[float, float, float], kind: str) -> tuple[float, float, float]:
    matrix = _CVD_MATRICES[kind]
    return tuple(  # type: ignore[return-value]
        min(1.0, max(0.0, sum(row[i] * linear[i] for i in range(3)))) for row in matrix
    )


def delta_e(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    a, b = oklab(left), oklab(right)
    return 100 * math.dist(a, b)


def relative_luminance(linear: tuple[float, float, float]) -> float:
    r, g, b = linear
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast_ratio(first: str, second: str) -> float:
    a, b = relative_luminance(linear_rgb(first)), relative_luminance(linear_rgb(second))
    lighter, darker = max(a, b), min(a, b)
    return (lighter + 0.05) / (darker + 0.05)


def validate(colours: list[str], mode: str, surface: str, all_pairs: bool) -> bool:
    low, high = LIGHTNESS_BAND[mode]
    ok = True
    print(f"palette of {len(colours)} colours, mode={mode}, surface={surface}\n")
    print(f"{'colour':>9}  {'L':>5} {'C':>5}  {'contrast':>8}  checks")
    for colour in colours:
        lab = oklab(linear_rgb(colour))
        lightness = lab[0]
        chroma = math.hypot(lab[1], lab[2])
        ratio = contrast_ratio(colour, surface)
        problems = []
        if not low <= lightness <= high:
            problems.append(f"lightness {lightness:.3f} outside {low}-{high}")
        if chroma < CHROMA_FLOOR:
            problems.append(f"chroma {chroma:.3f} < {CHROMA_FLOOR}")
        if ratio < CONTRAST_FLOOR:
            problems.append(f"contrast {ratio:.2f} < {CONTRAST_FLOOR} (needs labels)")
        ok &= not problems
        status = "PASS" if not problems else "FAIL: " + "; ".join(problems)
        print(f"{colour:>9}  {lightness:5.3f} {chroma:5.3f}  {ratio:8.2f}  {status}")

    pairs = (
        list(itertools.combinations(range(len(colours)), 2))
        if all_pairs
        else [(i, i + 1) for i in range(len(colours) - 1)]
    )
    label = "all pairs" if all_pairs else "adjacent pairs"
    print(f"\nseparation over {label} ({len(pairs)} pairs)")
    worst_normal = math.inf
    worst_cvd = math.inf
    for i, j in pairs:
        left, right = linear_rgb(colours[i]), linear_rgb(colours[j])
        normal = delta_e(left, right)
        worst_normal = min(worst_normal, normal)
        row = f"  {colours[i]} vs {colours[j]}  normal {normal:6.1f}"
        for kind in _CVD_MATRICES:
            value = delta_e(simulate(left, kind), simulate(right, kind))
            worst_cvd = min(worst_cvd, value)
            row += f"  {kind[:5]} {value:6.1f}"
        flags = []
        if normal < NORMAL_VISION_FLOOR:
            flags.append("NORMAL-VISION FAIL")
        pair_cvd = min(
            delta_e(simulate(left, kind), simulate(right, kind)) for kind in _CVD_MATRICES
        )
        if pair_cvd < CVD_FLOOR:
            flags.append("CVD FAIL")
        elif pair_cvd < CVD_TARGET:
            flags.append("cvd below target (needs secondary encoding)")
        if flags:
            ok = False if "FAIL" in " ".join(flags) else ok
            row += "   <- " + "; ".join(flags)
        print(row)
    print(
        f"\nworst normal-vision dE {worst_normal:.1f} (floor {NORMAL_VISION_FLOOR}); "
        f"worst CVD dE {worst_cvd:.1f} (target {CVD_TARGET}, floor {CVD_FLOOR})"
    )
    if worst_normal < NORMAL_VISION_FLOOR or worst_cvd < CVD_FLOOR:
        ok = False
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("palette", help="comma-separated hex colours")
    parser.add_argument("--mode", choices=("light", "dark"), default="light")
    parser.add_argument("--surface")
    parser.add_argument(
        "--pairs",
        choices=("adjacent", "all"),
        default="adjacent",
        help="'all' for scatter/small-multiples, where any two marks may touch",
    )
    args = parser.parse_args()
    colours = [item.strip() for item in args.palette.split(",") if item.strip()]
    surface = args.surface or DEFAULT_SURFACE[args.mode]
    return 0 if validate(colours, args.mode, surface, args.pairs == "all") else 1


if __name__ == "__main__":
    sys.exit(main())
