# tests/test_geno.py
import math
from pathlib import Path
import pytest

from src import Geno


def write_temp(dirpath: Path, name: str, content: str) -> Path:
    p = dirpath / name
    p.write_text((content.strip() + "\n"), encoding="utf-8")
    return p


# ---------- SEG × SEG ----------

def test_seg_seg_basic_overlap_example(tmp_path: Path):
    g = Geno()
    a = write_temp(tmp_path, "a.s", """
        1 2
        3 6
    """)
    b = write_temp(tmp_path, "b.s", """
        0 1
        1 5
    """)
    overlap = g.calBothSeg(str(a), str(b))
    assert overlap == 3  # positions 1, 3, 4


def test_seg_seg_touching_but_no_overlap(tmp_path: Path):
    g = Geno()
    a = write_temp(tmp_path, "a.s", "100 200")
    b = write_temp(tmp_path, "b.s", "200 300")
    overlap = g.calBothSeg(str(a), str(b))
    assert overlap == 0


def test_seg_seg_full_containment(tmp_path: Path):
    g = Geno()
    a = write_temp(tmp_path, "a.s", "10 50")
    b = write_temp(tmp_path, "b.s", "20 30")
    overlap = g.calBothSeg(str(a), str(b))
    assert overlap == 10


# ---------- FUNC × FUNC ----------

def test_func_func_correlation_matches_known_sample(tmp_path: Path):
    g = Geno()
    x = write_temp(tmp_path, "x.f", """
        10.0
        11.0
        12.0
        13.0
        14.0
        15.0
        16.0
    """)
    y = write_temp(tmp_path, "y.f", """
        10.5
        11.5
        12.0
        13.0
        13.5
        15.0
        14.0
    """)
    r = g.calBothFun(str(x), str(y))
    assert math.isclose(r, 0.9452853306994897, abs_tol=1e-9), f"r={r}"


def test_func_func_raises_when_one_series_is_constant_zero_std(tmp_path: Path):
    g = Geno()
    x = write_temp(tmp_path, "x.f", """
        1
        1
        1
        1
    """)
    y = write_temp(tmp_path, "y.f", """
        1
        2
        3
        4
    """)
    with pytest.raises(ValueError, match="Standard deviation is zero"):
        g.calBothFun(str(x), str(y))


def test_func_func_raises_when_lengths_differ(tmp_path: Path):
    g = Geno()
    short = write_temp(tmp_path, "x.f", """
        1
        2
        3
    """)
    long = write_temp(tmp_path, "y.f", """
        1
        2
        3
        4
    """)
    with pytest.raises(ValueError, match="Function files must have the same size"):
        g.calBothFun(str(short), str(long))


# ---------- SEG × FUNC ----------

def test_seg_func_mean_matches_known_sample(tmp_path: Path):
    g = Geno()
    seg = write_temp(tmp_path, "x.s", """
        1 2
        3 6
    """)
    fun_file = write_temp(tmp_path, "y.f", """
        10.5
        11.5
        12.0
        13.0
        13.5
        15.0
        14.0
    """)
    mean = g.calSegAndFun(str(seg), str(fun_file))
    assert math.isclose(mean, 13.25, abs_tol=1e-9), f"mean={mean}"


def test_seg_func_raises_when_no_covered_values(tmp_path: Path):
    g = Geno()
    seg = write_temp(tmp_path, "x.s", """
        10 11
        20 21
    """)
    fun_file = write_temp(tmp_path, "y.f", """
        1
        2
        3
    """)
    with pytest.raises(ValueError, match="No function values fall within the segments"):
        g.calSegAndFun(str(seg), str(fun_file))


# ---------- SEG file errors ----------

def test_seg_file_invalid_format_raises(tmp_path: Path):
    g = Geno()
    bad_seg = write_temp(tmp_path, "bad.s", """
        0 1 2
    """)
    other = write_temp(tmp_path, "ok.s", """
        5 6
    """)
    with pytest.raises(ValueError, match="Invalid line format"):
        g.calBothSeg(str(bad_seg), str(other))
