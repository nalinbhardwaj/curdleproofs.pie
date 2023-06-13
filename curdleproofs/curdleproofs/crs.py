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
)
from py_ecc.bls.g2_primitives import G1_to_pubkey, pubkey_to_G1
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
    def new(
        cls: Type[T_CurdleproofsCrs], ell: int, n_blinders: int
    ) -> T_CurdleproofsCrs:
        vec_G: List[PointProjective] = [get_random_point() for i in range(0, ell)]
        vec_H: List[PointProjective] = [
            get_random_point() for i in range(0, n_blinders)
        ]
        H: PointProjective = get_random_point()
        G_t: PointProjective = get_random_point()
        G_u: PointProjective = get_random_point()
        G_sum: PointProjective = reduce(add, vec_G, Z1)
        H_sum: PointProjective = reduce(add, vec_H, Z1)
        return cls(
            vec_G=vec_G,
            vec_H=vec_H,
            H=H,
            G_t=G_t,
            G_u=G_u,
            G_sum=G_sum,
            H_sum=H_sum,
        )

    def to_bytes(self) -> bytes:
        return b''.join([
            b''.join([G1_to_pubkey(g) for g in self.vec_G]),
            b''.join([G1_to_pubkey(h) for h in self.vec_H]),
            G1_to_pubkey(self.H),
            G1_to_pubkey(self.G_t),
            G1_to_pubkey(self.G_u),
            G1_to_pubkey(self.G_sum),
            G1_to_pubkey(self.H_sum),
        ])
    
    @classmethod
    def from_bytes(cls: Type[T_CurdleproofsCrs], b: bytes, ell: int, n_blinders: int) -> T_CurdleproofsCrs:
        def read_g1_point(i: int) -> PointProjective:
            return pubkey_to_G1(BLSPubkey(b[48 * i:48 * (i + 1)]))
        
        return cls(
            vec_G=[read_g1_point(i) for i in range(0, ell)],
            vec_H=[read_g1_point(ell + i) for i in range(0, n_blinders)],
            H=read_g1_point(ell + n_blinders),
            G_t=read_g1_point(ell + n_blinders + 1),
            G_u=read_g1_point(ell + n_blinders + 2),
            G_sum=read_g1_point(ell + n_blinders + 3),
            H_sum=read_g1_point(ell + n_blinders + 4),
        )

