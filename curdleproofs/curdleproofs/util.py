from random import randint
from math import log2
from typing import List, Tuple, Type, TypeVar, Union
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
from py_ecc.bls.hash import os2ip, i2osp
from py_ecc.bls.g2_primitives import G1_to_pubkey, pubkey_to_G1
from eth_typing import BLSPubkey
from py_ecc.bls.point_compression import compress_G1


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


def g1_from_bytes(b: bytes, offset_point: int) -> PointProjective:
    return pubkey_to_G1(BLSPubkey(b[48 * offset_point:48 * (offset_point + 1)]))


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


def field_to_json(f: FQ_type) -> str:
    return str(int(f))


T_JSON_FIELD = TypeVar("T_JSON_FIELD", Fr, FQ)


def field_from_json(s: str, field: Type[T_JSON_FIELD]) -> T_JSON_FIELD:
    return field(int(s))


def point_affine_to_json(p: PointAffine) -> Tuple[str, str]:
    return (field_to_json(p[0]), field_to_json(p[1]))


def point_projective_to_json(p: PointProjective) -> Tuple[str, str]:
    return point_affine_to_json(normalize(p))


def point_affine_from_json(t: Tuple[str, str]) -> PointAffine:
    return (field_from_json(t[0], FQ), field_from_json(t[1], FQ))


def point_projective_from_json(t: Tuple[str, str]) -> PointProjective:
    return affine_to_projective(point_affine_from_json(t))


def g1_to_bytes(p: PointProjective) -> bytes:
    return compress_G1(p).to_bytes(48, 'big')


def g1_list_to_bytes(ps: List[PointProjective]) -> bytes:
    return b''.join([g1_to_bytes(p) for p in ps])


def fr_to_bytes(fr: Fr) -> bytes:
    return fr.n.to_bytes(48, "big")


def log2_int(x: int) -> int:
    lg_x = int(log2(x))
    if x != 2**lg_x:
        raise Exception("x not a power of 2", x)
    return lg_x


class BufReader:
    def __init__(self, data):
        self.data = data
        self.ptr = 0

    def read_g1(self) -> PointProjective:
        end_ptr = self.ptr + 48
        p = pubkey_to_G1(BLSPubkey(self.data[self.ptr:end_ptr]))
        self.ptr = end_ptr
        return p
    
    def read_fr(self) -> Fr:
        end_ptr = self.ptr + 48
        p = Fr(os2ip(self.data[self.ptr:end_ptr]))
        self.ptr = end_ptr
        return p
