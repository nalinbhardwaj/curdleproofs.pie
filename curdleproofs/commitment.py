from audioop import mul
from typing import TypeVar
from curdleproofs.crs import CurdleproofsCrs
from util import PointAffine, PointProjective, Fr, field_to_bytes
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, Z1

T_GroupCommitment = TypeVar('T_GroupCommitment', bound="GroupCommitment")

class GroupCommitment:
  T_1: PointProjective
  T_2: PointProjective

  def __init__(self, crs_G: PointProjective, crs_H: PointProjective, T: PointProjective, r: Fr) -> None:
    self.T_1 = multiply(crs_G, r)
    self.T_2 = add(T, multiply(crs_H, r))
  
  def __init__(self, T_1: PointProjective, T_2: PointProjective) -> None:
    self.T_1 = T_1
    self.T_2 = T_2

  def __add__(self, other: T_GroupCommitment) -> T_GroupCommitment:
    return GroupCommitment(add(self.T_1, other.T_1), add(self.T_2, other.T_2))
  
  def __mul__(self, other: Fr) -> T_GroupCommitment:
    return GroupCommitment(multiply(self.T_1, other), multiply(self.T_2, other))
