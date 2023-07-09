from typing import Type, TypeVar
from curdleproofs.util import (
    point_projective_from_json,
    point_projective_to_json,
    BufReader,
    g1_to_bytes,
)
from py_arkworks_bls12381 import G1Point, Scalar


T_GroupCommitment = TypeVar("T_GroupCommitment", bound="GroupCommitment")


class GroupCommitment:
    T_1: G1Point
    T_2: G1Point

    def __init__(self, T_1: G1Point, T_2: G1Point) -> None:
        self.T_1 = T_1
        self.T_2 = T_2

    @classmethod
    def new(
        cls: Type[T_GroupCommitment],
        crs_G: G1Point,
        crs_H: G1Point,
        T: G1Point,
        r: Scalar,
    ) -> T_GroupCommitment:
        return cls(crs_G * r, T + crs_H * r)

    def __add__(self: T_GroupCommitment, other: object) -> T_GroupCommitment:
        if not isinstance(other, GroupCommitment):
            return NotImplemented
        return type(self)(self.T_1 + other.T_1, self.T_2 + other.T_2)

    def __mul__(self: T_GroupCommitment, other: object) -> T_GroupCommitment:
        if not isinstance(other, Scalar):
            return NotImplemented
        return type(self)(
            self.T_1 * other, self.T_2 * other
        )

    def __eq__(self: T_GroupCommitment, __o: object) -> bool:
        if not isinstance(__o, GroupCommitment):
            return NotImplemented
        return self.T_1 == __o.T_1 and self.T_2 == __o.T_2

    def to_json(self):
        return {
            "T_1": point_projective_to_json(self.T_1),
            "T_2": point_projective_to_json(self.T_2),
        }

    @classmethod
    def from_json(cls: Type[T_GroupCommitment], json) -> T_GroupCommitment:
        return cls(
            T_1=point_projective_from_json(json["T_1"]),
            T_2=point_projective_from_json(json["T_2"]),
        )

    def to_bytes(self) -> bytes:
        return b''.join([
            g1_to_bytes(self.T_1),
            g1_to_bytes(self.T_2),
        ])

    @classmethod
    def from_bytes(cls: Type[T_GroupCommitment], b: BufReader) -> T_GroupCommitment:
        return cls(
            T_1=b.read_g1(),
            T_2=b.read_g1(),
        )
