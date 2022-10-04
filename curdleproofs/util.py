from random import randint
from typing import List, TypeVar, Union
from py_ecc.typing import (
    Optimized_Field,
    Optimized_Point2D,
    Optimized_Point3D,
    FQ as FQ_type,
)
from py_ecc.optimized_bls12_381.optimized_curve import (
    curve_order,
    G1,
    multiply,
    normalize,
    add,
    neg,
    FQ,
)


class Fr(FQ_type):
    field_modulus: int = curve_order


PointAffine = Optimized_Point2D[Optimized_Field]
PointProjective = Optimized_Point3D[Optimized_Field]


def point_affine_to_bytes(point: PointAffine) -> bytes:
    return point[0].n.to_bytes(48, "big") + point[1].n.to_bytes(48, "big")


def points_affine_to_bytes(points: List[PointAffine]) -> List[bytes]:
    return [point_affine_to_bytes(point) for point in points]


def point_projective_to_bytes(point: PointProjective) -> bytes:
    return point_affine_to_bytes(normalize(point))


def points_projective_to_bytes(points: List[PointProjective]) -> List[bytes]:
    return [point_projective_to_bytes(point) for point in points]


def field_to_bytes(field: Fr) -> bytes:
    return field.n.to_bytes(48, "big")


def fields_to_bytes(fields: List[Fr]) -> List[bytes]:
    return [field_to_bytes(field) for field in fields]


def affine_to_projective(point: PointAffine) -> PointProjective:
    return (point[0], point[1], FQ.one())


def invert(f: Fr) -> Fr:
    res = Fr.one() / f
    assert res * f == Fr.one()  # fail in case f == 0
    return res


def get_random_point() -> PointProjective:
    a = randint(1, curve_order - 1)
    return multiply(G1, a)


def get_verification_scalars_bitstring(n: int, lg_n: int) -> List[List[int]]:
    bitstrings: List[List[int]] = []

    for i in range(0, n):
        bs = bin(i)[2:].zfill(lg_n)
        bitstrings.append([j for j in range(0, lg_n) if bs[j] == "1"])

    return bitstrings


def generate_blinders(n: int) -> List[Fr]:
    return [Fr(randint(0, Fr.field_modulus)) for _ in range(0, n)]


def inner_product(a: List[Fr], b: List[Fr]) -> Fr:
    assert len(a) == len(b)
    return sum([a[i] * b[i] for i in range(0, len(a))], Fr.zero())


T_GET_PERMUTATION = TypeVar("T_GET_PERMUTATION")


def get_permutation(
    vec_a: List[T_GET_PERMUTATION], permutation: Union[List[Fr], List[int]]
) -> List[T_GET_PERMUTATION]:
    return [vec_a[int(i)] for i in permutation]
