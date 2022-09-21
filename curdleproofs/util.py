from typing import List
from py_ecc.typing import (
    Optimized_Field,
    Optimized_Point2D,
    Optimized_Point3D,
    bls12_381_FQ,
)

Fr = bls12_381_FQ
PointAffine = Optimized_Point2D[Optimized_Field]
PointProjective = Optimized_Point3D[Optimized_Field]

def point_to_bytes(point: PointAffine) -> bytes:
    return point[0].to_bytes(48, 'big') + point[1].to_bytes(48, 'big')

def points_to_bytes(points: List[PointAffine]) -> bytes:
    return b''.join([point_to_bytes(point) for point in points])

def point_to_bytes(point: PointProjective) -> bytes:
    return point[0].to_bytes(48, 'big') + point[1].to_bytes(48, 'big') + point[2].to_bytes(48, 'big')

def points_to_bytes(points: List[PointProjective]) -> bytes:
    return b''.join([point_to_bytes(point) for point in points])

def field_to_bytes(field: Fr) -> bytes:
    return field.to_bytes(48, 'big')

def fields_to_bytes(fields: List[Fr]) -> bytes:
    return b''.join([field_to_bytes(field) for field in fields])