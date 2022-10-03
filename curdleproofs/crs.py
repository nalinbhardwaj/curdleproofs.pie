from functools import reduce
from random import randint
from typing import List
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, Z1
from util import PointAffine, PointProjective, affine_to_projective
import itertools

def get_random_point() -> PointProjective:
  a = randint(1, curve_order-1)
  return multiply(G1, a)

class CurdleproofsCrs:
  def __init__(self, ell: int, n_blinders: int) -> None:
    self.vec_G: List[PointAffine] = [normalize(get_random_point()) for i in range(0, ell)]
    self.vec_H: List[PointAffine] = [normalize(get_random_point()) for i in range(0, n_blinders)]
    self.H: PointProjective = get_random_point()
    self.G_t: PointProjective = get_random_point()
    self.G_u: PointProjective = get_random_point()
    self.G_sum: PointAffine = normalize(reduce(add, list(map(affine_to_projective, self.vec_G)), Z1))
    self.H_sum: PointAffine = normalize(reduce(add, list(map(affine_to_projective, self.vec_H)), Z1))

# crs = CurdleproofsCrs(3)
# print(crs.vec_G)
# print(crs.vec_H)
# print(crs.H)