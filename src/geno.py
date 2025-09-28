# geno.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List
import math


@dataclass(frozen=True, slots=True)
class _Seg:
    start: int
    end: int


class Geno:
    """
    Geno class provides methods to calculate overlaps and correlations
    between genomic segments (.s) and functions (.f).
    """

    # Reads a .s file and returns a list of _Seg objects.
    @staticmethod
    def readSegFile(filePath: str) -> List[_Seg]:
        path = Path(filePath)
        lines = path.read_text(encoding="utf-8").splitlines()

        segs: List[_Seg] = []
        for line in lines:
            parts = line.split()  # split on any whitespace
            if len(parts) != 2:
                raise ValueError(f"Invalid line format: {line}")
            try:
                start = int(parts[0])
                end = int(parts[1])
            except ValueError:
                raise ValueError(f"Invalid start/end position: {line}")

            if not (start < end):
                raise ValueError(f"Start position must be less than end position: {line}")

            segs.append(_Seg(start, end))
        return segs

    # Reads a .f file and returns a list of float values.
    @staticmethod
    def readFunFile(filePath: str) -> List[float]:
        path = Path(filePath)
        lines = path.read_text(encoding="utf-8").splitlines()

        values: List[float] = []
        for line in lines:
            try:
                values.append(float(line.strip()))
            except ValueError:
                raise ValueError(f"Invalid value: {line}")
        return values

    # Calculates the total overlap length between two segment files.
    def calBothSeg(self, file1: str, file2: str) -> int:
        segs1 = self.readSegFile(file1)
        segs2 = self.readSegFile(file2)

        overlapCnt = 0
        i = 0
        j = 0

        while i < len(segs1) and j < len(segs2):
            s = max(segs1[i].start, segs2[j].start)
            e = min(segs1[i].end, segs2[j].end)

            if s < e:
                overlapCnt += (e - s)

            if segs1[i].end < segs2[j].end:
                i += 1
            else:
                j += 1

        return overlapCnt

    # Calculates the Pearson correlation coefficient between two function files.
    def calBothFun(self, file1: str, file2: str) -> float:
        funs1 = self.readFunFile(file1)
        funs2 = self.readFunFile(file2)

        if len(funs1) != len(funs2):
            raise ValueError("Function files must have the same size.")

        # mean
        mean1 = sum(funs1) / len(funs1)
        mean2 = sum(funs2) / len(funs2)

        # sqrt of sum of squared deviations
        sqrtSquareDiff1 = math.sqrt(sum((x - mean1) * (x - mean1) for x in funs1))
        sqrtSquareDiff2 = math.sqrt(sum((y - mean2) * (y - mean2) for y in funs2))

        if sqrtSquareDiff1 == 0.0 or sqrtSquareDiff2 == 0.0:
            raise ValueError("Standard deviation is zero, cannot calculate correlation.")

        # numerator
        total = 0.0
        for x, y in zip(funs1, funs2):
            total += (x - mean1) * (y - mean2)

        return total / (sqrtSquareDiff1 * sqrtSquareDiff2)

    # Calculates the average function value within the segments defined in the segment file.
    def calSegAndFun(self, file1: str, file2: str) -> float:
        segs = self.readSegFile(file1)
        funs = self.readFunFile(file2)

        sum_vals = 0.0
        count = 0
        segIndex = 0

        for i in range(len(funs)):
            while segIndex < len(segs) and i >= segs[segIndex].end:
                segIndex += 1

            if segIndex < len(segs) and (segs[segIndex].start <= i < segs[segIndex].end):
                sum_vals += funs[i]
                count += 1

        if count == 0:
            raise ValueError("No function values fall within the segments.")

        return sum_vals / count


if __name__ == "__main__":
    print(Geno().calBothSeg("./data/x.s", "./data/y.s"))
    print(Geno().calBothSeg("./data/testfile_a.s", "./data/testfile_b.s"))

    print(Geno().calBothFun("./data/x.f", "./data/y.f"))
    print(Geno().calBothFun("./data/testfile_a.f", "./data/testfile_b.f"))

    print(Geno().calSegAndFun("./data/x.s", "./data/y.f"))
    print(Geno().calSegAndFun("./data/testfile_a.s", "./data/testfile_b.f"))
