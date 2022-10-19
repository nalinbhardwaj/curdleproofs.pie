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
