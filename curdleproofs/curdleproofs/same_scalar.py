import random
from curdleproofs.commitment import GroupCommitment
from curdleproofs.util import field_from_json, field_to_json, points_projective_to_bytes
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import Type, TypeVar
from curdleproofs.util import (
    PointProjective,
    Fr,
    BufReader,
    fr_to_bytes,
)
from py_ecc.optimized_bls12_381.optimized_curve import multiply

T_SameScalarProof = TypeVar("T_SameScalarProof", bound="SameScalarProof")


class SameScalarProof:
    def __init__(
        self, cm_A: GroupCommitment, cm_B: GroupCommitment, z_k: Fr, z_t: Fr, z_u: Fr
    ) -> None:
        self.cm_A = cm_A
        self.cm_B = cm_B
        self.z_k = z_k
        self.z_t = z_t
        self.z_u = z_u

    @classmethod
    def new(
        cls: Type[T_SameScalarProof],
        crs_G_t: PointProjective,
        crs_G_u: PointProjective,
        crs_H: PointProjective,
        R: PointProjective,
        S: PointProjective,
        cm_T: GroupCommitment,
        cm_U: GroupCommitment,
        k: Fr,
        r_t: Fr,
        r_u: Fr,
        transcript: CurdleproofsTranscript,
    ) -> T_SameScalarProof:
        r_a = Fr(random.randint(1, Fr.field_modulus))
        r_b = Fr(random.randint(1, Fr.field_modulus))
        r_k = Fr(random.randint(1, Fr.field_modulus))

        cm_A = GroupCommitment.new(crs_G_t, crs_H, multiply(R, int(r_k)), r_a)
        cm_B = GroupCommitment.new(crs_G_u, crs_H, multiply(S, int(r_k)), r_b)

        transcript.append_list(
            b"sameexp_points",
            points_projective_to_bytes(
                [
                    R,
                    S,
                    cm_T.T_1,
                    cm_T.T_2,
                    cm_U.T_1,
                    cm_U.T_2,
                    cm_A.T_1,
                    cm_A.T_2,
                    cm_B.T_1,
                    cm_B.T_2,
                ]
            ),
        )
        alpha = transcript.get_and_append_challenge(b"same_scalar_alpha")

        z_k = r_k + k * alpha
        z_t = r_a + r_t * alpha
        z_u = r_b + r_u * alpha

        return cls(cm_A, cm_B, z_k, z_t, z_u)

    def verify(
        self,
        crs_G_t: PointProjective,
        crs_G_u: PointProjective,
        crs_H: PointProjective,
        R: PointProjective,
        S: PointProjective,
        cm_T: GroupCommitment,
        cm_U: GroupCommitment,
        transcript: CurdleproofsTranscript,
    ):
        transcript.append_list(
            b"sameexp_points",
            points_projective_to_bytes(
                [
                    R,
                    S,
                    cm_T.T_1,
                    cm_T.T_2,
                    cm_U.T_1,
                    cm_U.T_2,
                    self.cm_A.T_1,
                    self.cm_A.T_2,
                    self.cm_B.T_1,
                    self.cm_B.T_2,
                ]
            ),
        )

        alpha = transcript.get_and_append_challenge(b"same_scalar_alpha")
        expected_1 = GroupCommitment.new(
            crs_G_t, crs_H, multiply(R, int(self.z_k)), self.z_t
        )
        expected_2 = GroupCommitment.new(
            crs_G_u, crs_H, multiply(S, int(self.z_k)), self.z_u
        )

        computed_1 = self.cm_A + (cm_T * alpha)
        computed_2 = self.cm_B + (cm_U * alpha)

        assert expected_1 == computed_1 and expected_2 == computed_2

    def to_json(self):
        return {
            "cm_A": self.cm_A.to_json(),
            "cm_B": self.cm_B.to_json(),
            "z_k": field_to_json(self.z_k),
            "z_t": field_to_json(self.z_t),
            "z_u": field_to_json(self.z_u),
        }

    @classmethod
    def from_json(cls: Type[T_SameScalarProof], json) -> T_SameScalarProof:
        return cls(
            cm_A=GroupCommitment.from_json(json["cm_A"]),
            cm_B=GroupCommitment.from_json(json["cm_B"]),
            z_k=field_from_json(json["z_k"], Fr),
            z_t=field_from_json(json["z_t"], Fr),
            z_u=field_from_json(json["z_u"], Fr),
        )

    def to_bytes(self) -> bytes:
        return b''.join([
            self.cm_A.to_bytes(),
            self.cm_B.to_bytes(),
            fr_to_bytes(self.z_k),
            fr_to_bytes(self.z_t),
            fr_to_bytes(self.z_u),
        ])

    @classmethod
    def from_bytes(cls: Type[T_SameScalarProof], b: BufReader) -> T_SameScalarProof:
        return cls(
            cm_A=GroupCommitment.from_bytes(b),
            cm_B=GroupCommitment.from_bytes(b),
            z_k=b.read_fr(),
            z_t=b.read_fr(),
            z_u=b.read_fr(),
        )
