import random
from curdleproofs.util import PointProjective, Fr, affine_to_projective
from typing import Dict, List, Tuple, Union
from py_ecc.optimized_bls12_381.optimized_curve import (
    curve_order,
    G1,
    multiply,
    normalize,
    add,
    Z1,
    eq,
    FQ,
    is_inf,
)


def compute_MSM(
    bases: List[PointProjective], scalars: Union[List[Fr], List[int]]
) -> PointProjective:
    current = Z1  # zero
    for (base, scalar) in zip(bases, scalars):
        current = add(current, multiply(base, int(scalar)))  # type: ignore
    return current


# class MSMAccumulatorInefficient:
#   # TODO: does not accumulate right now
#   def __init__(self) -> None:
#     self.MSMs: List[Tuple[SingleMSM, PointProjective]] = []

#   def accumulate_check(self, C: PointProjective, bases: List[PointProjective], scalars: List[int]) -> None:
#     # print("accumulating", C, bases, scalars)
#     self.MSMs.append((SingleMSM(bases, scalars), C))

#   def verify(self):
#     for msm in self.MSMs:
#       # print("computed", normalize(msm[0].compute()), "expected", normalize(msm[1]), "eq", eq(msm[0].compute(), msm[1]))
#       if not eq(msm[0].compute(), msm[1]):
#         return False
#     return True


class MSMAccumulator:
    def __init__(self) -> None:
        self.A_c = Z1
        self.base_scalar_map: Dict[Tuple[int, int], Fr] = {}

    def accumulate_check(
        self,
        C: PointProjective,
        bases: List[PointProjective],
        scalars: Union[List[Fr], List[int]],
    ) -> None:
        random_factor = Fr(random.randint(1, Fr.field_modulus))

        self.A_c = add(self.A_c, multiply(C, int(random_factor)))

        for (base, scalar) in zip(bases, scalars):
            # print("base", base)
            if is_inf(base):
                continue
            base_affine_int_untyped = tuple(map(int, normalize(base)))
            base_affine_int = (base_affine_int_untyped[0], base_affine_int_untyped[1])
            # print("base_affine_int", base_affine_int)
            if base_affine_int not in self.base_scalar_map:
                self.base_scalar_map[base_affine_int] = Fr.zero()
            self.base_scalar_map[base_affine_int] = self.base_scalar_map[
                base_affine_int
            ] + random_factor * Fr(scalar)

    def verify(self):
        bases: List[Tuple[int, int]]
        scalars: List[Fr]
        bases, scalars = map(list, zip(*self.base_scalar_map.items()))  # type: ignore
        computed = compute_MSM(
            list(map(lambda t: affine_to_projective((FQ(t[0]), FQ(t[1]))), bases)),
            list(map(int, scalars)),
        )
        # print("bases", bases, "scalars", scalars, "computed", normalize(computed), "expected", normalize(self.A_c), "eq", eq(computed, self.A_c))
        assert eq(computed, self.A_c)
