from util import PointAffine, PointProjective
from typing import List
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, Z1


class SingleMSM:
  def __init__(self, C: PointProjective, bases: List[PointProjective], scalars: List[int]) -> None:
    self.C = C
    self.bases = bases
    self.scalars = scalars

  def check(self) -> bool:
    current = Z1 # zero
    for (base, scalar) in zip(self.bases, self.scalars):
      current = add(current, multiply(base, scalar))
    return self.C == current



class MSMAccumulatorInefficient:
  # TODO: does not accumulate right now
  def __init__(self) -> None:
    self.MSMs: List[SingleMSM] = []
  
  def accumulate_check(self, C: PointProjective, bases: List[PointProjective], scalars: List[int]) -> None:
    self.MSMs.append(SingleMSM(C, bases, scalars))

  def verify(self):
    for msm in self.MSMs:
      if not msm.check():
        return False
    return True