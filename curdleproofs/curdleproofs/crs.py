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
        count = ell + n_blinders + 3
        points: List[PointProjective] = [get_random_point() for _ in range(0, count)]
        return cls.from_random_points(ell, n_blinders, points)

    @classmethod
    def from_random_points(cls: Type[T_CurdleproofsCrs], ell: int, n_blinders: int, points: List[PointProjective]) -> T_CurdleproofsCrs:
        vec_G = points[0:ell]
        vec_H = points[ell:ell+n_blinders]
        return cls(
            vec_G=vec_G,
            vec_H=vec_H,
            H=points[ell+n_blinders + 0],
            G_t=points[ell+n_blinders + 1],
            G_u=points[ell+n_blinders + 2],
            G_sum=reduce(add, vec_G, Z1),
            H_sum=reduce(add, vec_H, Z1),
        )

    def to_json(self) -> str:
        dic = {
            "vec_G": [point_projective_to_json(g) for g in self.vec_G],
            "vec_H": [point_projective_to_json(h) for h in self.vec_H],
            "H": point_projective_to_json(self.H),
            "G_t": point_projective_to_json(self.G_t),
            "G_u": point_projective_to_json(self.G_u),
            "G_sum": point_projective_to_json(self.G_sum),
            "H_sum": point_projective_to_json(self.H_sum),
        }
        return json.dumps(dic)

    @classmethod
    def from_json(cls: Type[T_CurdleproofsCrs], json_str: str) -> T_CurdleproofsCrs:
        dic = json.loads(json_str)
        return cls(
            vec_G=[point_projective_from_json(g) for g in dic["vec_G"]],
            vec_H=[point_projective_from_json(h) for h in dic["vec_H"]],
            H=point_projective_from_json(dic["H"]),
            G_t=point_projective_from_json(dic["G_t"]),
            G_u=point_projective_from_json(dic["G_u"]),
            G_sum=point_projective_from_json(dic["G_sum"]),
            H_sum=point_projective_from_json(dic["H_sum"]),
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
    def from_bytes(cls: Type[T_CurdleproofsCrs], b: BufReader, ell: int, n_blinders: int) -> T_CurdleproofsCrs:
        return cls(
            vec_G=[b.read_g1() for _ in range(0, ell)],
            vec_H=[b.read_g1() for _ in range(0, n_blinders)],
            H=b.read_g1(),
            G_t=b.read_g1(),
            G_u=b.read_g1(),
            G_sum=b.read_g1(),
            H_sum=b.read_g1(),
        )
