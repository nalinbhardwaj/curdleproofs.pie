import json
import random
from typing import Any, Dict, Tuple, Type, TypeVar
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.util import (
    PointAffine,
    PointProjective,
    Fr,
    field_to_bytes,
    get_random_point,
    point_projective_from_json,
    point_projective_to_json,
)
from py_ecc.optimized_bls12_381.optimized_curve import (
    curve_order,
    G1,
    multiply,
    normalize,
    add,
    Z1,
    eq,
)

T_GroupCommitment = TypeVar("T_GroupCommitment", bound="GroupCommitment")


class GroupCommitment:
    T_1: PointProjective
    T_2: PointProjective

    def __init__(self, T_1: PointProjective, T_2: PointProjective) -> None:
        self.T_1 = T_1
        self.T_2 = T_2

    @classmethod
    def new(
        cls: Type[T_GroupCommitment],
        crs_G: PointProjective,
        crs_H: PointProjective,
        T: PointProjective,
        r: Fr,
    ) -> T_GroupCommitment:
        return cls(multiply(crs_G, int(r)), add(T, multiply(crs_H, int(r))))

    def __add__(self: T_GroupCommitment, other: object) -> T_GroupCommitment:
        if not isinstance(other, GroupCommitment):
            return NotImplemented
        return type(self)(add(self.T_1, other.T_1), add(self.T_2, other.T_2))

    def __mul__(self: T_GroupCommitment, other: object) -> T_GroupCommitment:
        if not isinstance(other, Fr):
            return NotImplemented
        return type(self)(
            multiply(self.T_1, int(other)), multiply(self.T_2, int(other))
        )

    def __eq__(self: T_GroupCommitment, __o: object) -> bool:
        if not isinstance(__o, GroupCommitment):
            return NotImplemented
        return eq(self.T_1, __o.T_1) and eq(self.T_2, __o.T_2)

    def to_json(self) -> str:
        dic = {
            "T_1": point_projective_to_json(self.T_1),
            "T_2": point_projective_to_json(self.T_2),
        }
        return json.dumps(dic)

    @classmethod
    def from_json(cls: Type[T_GroupCommitment], json_str: str) -> T_GroupCommitment:
        dic = json.loads(json_str)
        return cls(
            T_1=point_projective_from_json(dic["T_1"]),
            T_2=point_projective_from_json(dic["T_2"]),
        )
