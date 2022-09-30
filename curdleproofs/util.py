from typing import List
from py_ecc.typing import (
    Optimized_Field,
    Optimized_Point2D,
    Optimized_Point3D,
    FQ,
)
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, neg

class Fr(FQ):
    field_modulus: int = curve_order

PointAffine = Optimized_Point2D[Optimized_Field]
PointProjective = Optimized_Point3D[Optimized_Field]

def point_affine_to_bytes(point: PointAffine) -> bytes:
    return point[0].n.to_bytes(48, 'big') + point[1].n.to_bytes(48, 'big')

def points_affine_to_bytes(points: List[PointAffine]) -> List[bytes]:
    return [point_affine_to_bytes(point) for point in points]

def point_projective_to_bytes(point: PointProjective) -> bytes:
    return point_affine_to_bytes(normalize(point))

def points_projective_to_bytes(points: List[PointProjective]) -> List[bytes]:
    return [point_projective_to_bytes(point) for point in points]

def field_to_bytes(field: Fr) -> bytes:
    return field.n.to_bytes(48, 'big')

def fields_to_bytes(fields: List[Fr]) -> List[bytes]:
    return [field_to_bytes(field) for field in fields]

def affine_to_projective(point: PointAffine) -> PointProjective:
    return (point[0], point[1], 1)

def invert(f: Fr) -> Fr:
    res = Fr.one() / f
    assert res * f == Fr.one() # fail in case f == 0
    return res