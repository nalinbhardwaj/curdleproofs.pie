from random import randint
from math import log2
from typing import List, TypeVar, NewType
from py_arkworks_bls12381 import G1Point, Scalar


CURVE_ORDER = 52435875175126190479447740508185965837690552500527637822603658699938581184513
# Generator
G1 = G1Point()
# Point at infinity over FQ
Z1 = G1Point.identity()


BLSPubkey = NewType('BLSPubkey', bytes)  # bytes48


def g1_is_inf(point: G1Point):
    return point == Z1


def random_scalar() -> Scalar:
    # Note: the constructor 'Scalar()' errors with integers of more than 128 bits
    # Scalar.from_le_bytes() requires integers less than  CURVE_ORDER
    return Scalar.from_le_bytes(randint(1, CURVE_ORDER - 1).to_bytes(32, 'little'))


def point_projective_to_bytes(point: G1Point) -> bytes:
    return bytes(point.to_compressed_bytes())


def points_projective_to_bytes(points: List[G1Point]) -> List[bytes]:
    return [point_projective_to_bytes(point) for point in points]


def point_projective_from_bytes(b: bytes) -> G1Point:
    return G1Point.from_compressed_bytes_unchecked(b)


def field_to_bytes(field: Scalar) -> bytes:
    return bytes(field.to_le_bytes())


def fields_to_bytes(fields: List[Scalar]) -> List[bytes]:
    return [field_to_bytes(field) for field in fields]


def g1_from_bytes(b: bytes, offset_point: int) -> G1Point:
    return point_projective_from_bytes(BLSPubkey(b[48 * offset_point:48 * (offset_point + 1)]))


def invert(f: Scalar) -> Scalar:
    res = f.inverse()
    assert res * f == Scalar(1)  # fail in case f == 0
    return res


def scalar_pow(f: Scalar, n: int) -> Scalar:
    result = Scalar(1)
    while n != 0:
        if n % 2 == 1:
            result *= f
        f *= f
        n //= 2
    return result


def get_random_point() -> G1Point:
    return G1 * random_scalar()


def get_verification_scalars_bitstring(n: int, lg_n: int) -> List[List[int]]:
    bitstrings: List[List[int]] = []

    for i in range(0, n):
        bs = bin(i)[2:].zfill(lg_n)
        bitstrings.append([j for j in range(0, lg_n) if bs[j] == "1"])

    return bitstrings


def generate_blinders(n: int) -> List[Scalar]:
    return [random_scalar() for _ in range(0, n)]


def inner_product(a: List[Scalar], b: List[Scalar]) -> Scalar:
    assert len(a) == len(b)
    return sum([a[i] * b[i] for i in range(0, len(a))], Scalar(0))


T_GET_PERMUTATION = TypeVar("T_GET_PERMUTATION")


def get_permutation(
    vec_a: List[T_GET_PERMUTATION], permutation: List[int]
) -> List[T_GET_PERMUTATION]:
    return [vec_a[int(i)] for i in permutation]


def field_to_json(f: Scalar) -> str:
    return bytes(f.to_le_bytes()).hex()


def field_from_json(s: str) -> Scalar:
    return Scalar.from_le_bytes(bytes.fromhex(s))


def scalar_from_bytes(b: bytes) -> Scalar:
    return Scalar.from_le_bytes(b)


def point_projective_to_json(p: G1Point) -> str:
    return bytes(p.to_compressed_bytes()).hex()


def point_projective_from_json(p_hex: str) -> G1Point:
    return G1Point.from_compressed_bytes_unchecked(bytes.fromhex(p_hex))


def g1_to_bytes(p: G1Point) -> bytes:
    return bytes(p.to_compressed_bytes())


def g1_list_to_bytes(ps: List[G1Point]) -> bytes:
    return b''.join([g1_to_bytes(p) for p in ps])


def fr_to_bytes(fr: Scalar) -> bytes:
    return field_to_bytes(fr)


def log2_int(x: int) -> int:
    lg_x = int(log2(x))
    if x != 2**lg_x:
        raise Exception("x not a power of 2", x)
    return lg_x


class BufReader:
    def __init__(self, data):
        self.data = data
        self.ptr = 0

    def read_g1(self) -> G1Point:
        end_ptr = self.ptr + 48
        p = point_projective_from_bytes(BLSPubkey(self.data[self.ptr:end_ptr]))
        self.ptr = end_ptr
        return p

    def read_fr(self) -> Scalar:
        end_ptr = self.ptr + 32
        p = scalar_from_bytes(self.data[self.ptr:end_ptr])
        self.ptr = end_ptr
        return p
