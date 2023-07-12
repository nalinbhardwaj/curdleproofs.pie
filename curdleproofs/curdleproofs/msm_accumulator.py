from curdleproofs.util import random_scalar, g1_is_inf, Z1, point_projective_to_bytes
from typing import Dict, List, Tuple
from py_arkworks_bls12381 import G1Point, Scalar


def compute_MSM(
    bases: List[G1Point], scalars: List[Scalar]
) -> G1Point:
    current = G1Point.identity()  # zero
    for (base, scalar) in zip(bases, scalars):
        current = current + (base * scalar)  # type: ignore
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
        self.base_scalar_map: Dict[bytes, Scalar] = {}

    def accumulate_check(
        self,
        C: G1Point,
        bases: List[G1Point],
        scalars: List[Scalar],
    ) -> None:
        random_factor = random_scalar()

        self.A_c = self.A_c + C * random_factor

        for (base, scalar) in zip(bases, scalars):
            # note: optimization, zero bases contribute nothing to the MSM
            if g1_is_inf(base):
                continue

            # Note: G1Point is not hashable so a different representation is necessary to index base_scalar_map
            # TODO: Compressing and decompressing each point is unnecessary, find a different hashble representation
            base_comp = point_projective_to_bytes(base)
            # print("base_affine_int", base_affine_int)
            if base_comp not in self.base_scalar_map:
                self.base_scalar_map[base_comp] = Scalar(0)
            self.base_scalar_map[base_comp] = self.base_scalar_map[base_comp] + random_factor * scalar

    def verify(self):
        bases: List[Tuple[int, int]]
        scalars: List[Scalar]
        bases, scalars = map(list, zip(*self.base_scalar_map.items()))  # type: ignore
        computed = compute_MSM(
            list(map(lambda t: G1Point.from_compressed_bytes_unchecked(t), bases)),
            scalars,
        )
        assert computed == self.A_c
