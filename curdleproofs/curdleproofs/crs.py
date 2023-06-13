from functools import reduce
import json
from random import randint
from typing import List, Type, TypeVar
from py_ecc.optimized_bls12_381.optimized_curve import (
    curve_order,
    G1,
    multiply,
    normalize,
    add,
    Z1,
)
from curdleproofs.util import (
    PointProjective,
    affine_to_projective,
    get_random_point,
    point_projective_from_json,
    point_projective_to_json,
    BufReader,
    g1_to_bytes,
)
from py_ecc.bls.g2_primitives import G1_to_pubkey, pubkey_to_G1
from py_ecc.bls.point_compression import decompress_G1
from py_ecc.bls.typing import G1Compressed
from eth_typing import BLSPubkey

T_CurdleproofsCrs = TypeVar("T_CurdleproofsCrs", bound="CurdleproofsCrs")


class CurdleproofsCrs:
    def __init__(
        self,
        vec_G: List[PointProjective],
        vec_H: List[PointProjective],
        H: PointProjective,
        G_t: PointProjective,
        G_u: PointProjective,
        G_sum: PointProjective,
        H_sum: PointProjective,
    ) -> None:
        self.vec_G = vec_G
        self.vec_H = vec_H
        self.H = H
        self.G_t = G_t
        self.G_u = G_u
        self.G_sum = G_sum
        self.H_sum = H_sum

    @classmethod
    def new(cls: Type[T_CurdleproofsCrs], n: int, n_blinders: int) -> T_CurdleproofsCrs:
        count = min_poins_required(n)
        points: List[PointProjective] = [get_random_point() for _ in range(0, count)]
        return cls.from_points(n, n_blinders, points)
    
    @classmethod
    def from_points_compressed(cls: Type[T_CurdleproofsCrs], n: int, n_blinders: int, points_compressed: List[BLSPubkey]) -> T_CurdleproofsCrs:
        count = min_poins_required(n)
        points = [pubkey_to_G1(BLSPubkey(p)) for p in points_compressed[0:min(count, len(points_compressed))]]
        return cls.from_points(n, n_blinders, points)
    
    @classmethod
    def from_points(cls: Type[T_CurdleproofsCrs], n: int, n_blinders: int, points: List[PointProjective]) -> T_CurdleproofsCrs:
        count = min_poins_required(n)
        if len(points) < count:
            raise Exception("not min points required", count, len(points))
        if n_blinders >= n:
            raise Exception("n_blinders >= n")

        vec_G = points[0:n - n_blinders]
        vec_H = points[n - n_blinders:n]
        return cls(
            vec_G=vec_G,
            vec_H=vec_H,
            H=points[n + 0],
            G_t=points[n + 1],
            G_u=points[n + 2],
            G_sum=reduce(add, vec_G, Z1),
            H_sum=reduce(add, vec_H, Z1),
        )
    
    def to_bytes(self) -> bytes:
        return b''.join([
            b''.join([g1_to_bytes(g) for g in self.vec_G]),
            b''.join([g1_to_bytes(h) for h in self.vec_H]),
            g1_to_bytes(self.H),
            g1_to_bytes(self.G_t),
            g1_to_bytes(self.G_u),
            g1_to_bytes(self.G_sum),
            g1_to_bytes(self.H_sum),
        ])
    
    @classmethod
    def from_bytes(cls: Type[T_CurdleproofsCrs], b: BufReader, n: int, n_blinders: int) -> T_CurdleproofsCrs:
        return cls(
            vec_G=[b.read_g1() for _ in range(0, n - n_blinders)],
            vec_H=[b.read_g1() for _ in range(0, n_blinders)],
            H=b.read_g1(),
            G_t=b.read_g1(),
            G_u=b.read_g1(),
            G_sum=b.read_g1(),
            H_sum=b.read_g1(),
        )

def min_poins_required(n: int) -> int:
    return n + 3