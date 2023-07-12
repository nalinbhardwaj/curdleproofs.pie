from functools import reduce
import json
from math import log2
from typing import List, Type, TypeVar
from curdleproofs.util import (
    get_random_point,
    point_projective_from_json,
    point_projective_to_json,
    BufReader,
    g1_to_bytes,
    Z1
)
from py_arkworks_bls12381 import G1Point


T_CurdleproofsCrs = TypeVar("T_CurdleproofsCrs", bound="CurdleproofsCrs")


class CurdleproofsCrs:
    def __init__(
        self,
        vec_G: List[G1Point],
        vec_H: List[G1Point],
        H: G1Point,
        G_t: G1Point,
        G_u: G1Point,
        G_sum: G1Point,
        H_sum: G1Point,
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
        points: List[G1Point] = [get_random_point() for _ in range(0, count)]
        return cls.from_random_points(ell, n_blinders, points)

    @classmethod
    def from_random_points(cls: Type[T_CurdleproofsCrs], ell: int, n_blinders: int, points: List[G1Point]) -> T_CurdleproofsCrs:
        min_points_count = ell + n_blinders + 3
        if len(points) < min_points_count:
            raise Exception("not min points required", min_points_count, len(points))

        n = ell + n_blinders
        if n != 2**int(log2(n)):
            raise Exception("ell + n_blinders not a power of 2, ell={} n_blinders={}".format(ell, n_blinders))

        vec_G = points[0:ell]
        vec_H = points[ell:ell + n_blinders]
        return cls(
            vec_G=vec_G,
            vec_H=vec_H,
            H=points[ell + n_blinders + 0],
            G_t=points[ell + n_blinders + 1],
            G_u=points[ell + n_blinders + 2],
            G_sum=reduce(lambda a, b: a + b, vec_G, Z1),
            H_sum=reduce(lambda a, b: a + b, vec_H, Z1),
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
