import random
from typing import TypeVar
from crs import CurdleproofsCrs, get_random_point
from util import PointAffine, PointProjective, Fr, field_to_bytes
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, Z1, eq

T_GroupCommitment = TypeVar('T_GroupCommitment', bound="GroupCommitment")

class GroupCommitment:
  T_1: PointProjective
  T_2: PointProjective

  def __init__(self, T_1: PointProjective, T_2: PointProjective) -> None:
    self.T_1 = T_1
    self.T_2 = T_2

  @classmethod
  def new(cls,
    crs_G: PointProjective,
    crs_H: PointProjective,
    T: PointProjective,
    r: Fr
  ) -> T_GroupCommitment:
    return cls(multiply(crs_G, int(r)), add(T, multiply(crs_H, int(r))))

  def __add__(self, other: object) -> T_GroupCommitment:
    if not isinstance(other, GroupCommitment):
      return NotImplemented
    return type(self)(add(self.T_1, other.T_1), add(self.T_2, other.T_2))
  
  def __mul__(self, other: object) -> T_GroupCommitment:
    if not isinstance(other, Fr):
      return NotImplemented
    return type(self)(multiply(self.T_1, int(other)), multiply(self.T_2, int(other)))

  def __eq__(self, __o: object) -> bool:
    if not isinstance(__o, GroupCommitment):
      return NotImplemented
    return eq(self.T_1, __o.T_1) and eq(self.T_2, __o.T_2)


def test_group_commit():
  crs_G = get_random_point()
  crs_H = get_random_point()

  A = get_random_point()
  B = get_random_point()

  r_a = Fr(random.randint(1, Fr.field_modulus))
  r_b = Fr(random.randint(1, Fr.field_modulus))

  cm_a = GroupCommitment.new(crs_G, crs_H, A, r_a)
  cm_b = GroupCommitment.new(crs_G, crs_H, B, r_b)
  cm_a_b = GroupCommitment.new(crs_G, crs_H, add(A, B), r_a + r_b)

  assert cm_a + cm_b == cm_a_b

# test_group_commit()