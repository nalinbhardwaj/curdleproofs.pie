from util import PointAffine, PointProjective
from typing import List, Tuple
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, Z1, eq


class SingleMSM:
  def __init__(self, bases: List[PointProjective], scalars: List[int]) -> None:
    self.bases = bases
    self.scalars = scalars
  
  def compute(self) -> PointProjective:
    current = Z1 # zero
    for (base, scalar) in zip(self.bases, self.scalars):
      current = add(current, multiply(base, scalar))
    return current


class MSMAccumulatorInefficient:
  # TODO: does not accumulate right now
  def __init__(self) -> None:
    self.MSMs: List[Tuple[SingleMSM, PointProjective]] = []
  
  def accumulate_check(self, C: PointProjective, bases: List[PointProjective], scalars: List[int]) -> None:
    # print("accumulating", C, bases, scalars)
    self.MSMs.append((SingleMSM(bases, scalars), C))

  def verify(self):
    for msm in self.MSMs:
      # print("computed", normalize(msm[0].compute()), "expected", normalize(msm[1]), "eq", eq(msm[0].compute(), msm[1]))
      if not eq(msm[0].compute(), msm[1]):
        return False
    return True