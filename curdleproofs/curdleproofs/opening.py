from typing import Type, TypeVar
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript

from curdleproofs.util import (
    field_from_json,
    field_to_json,
    generate_blinders,
    point_projective_from_json,
    point_projective_to_json,
    points_projective_to_bytes,
    BufReader,
    g1_to_bytes,
    fr_to_bytes,
    G1,
)
from py_arkworks_bls12381 import G1Point, Scalar

T_TrackerOpeningProof = TypeVar("T_TrackerOpeningProof", bound="TrackerOpeningProof")


# proof of knowledge of `k` s.t. k_r_G = k*r_G and k_G = k*G
class TrackerOpeningProof:
    def __init__(
        self,
        A: G1Point,
        B: G1Point,
        s: Scalar,
    ) -> None:
        self.A = A
        self.B = B
        self.s = s

    @classmethod
    def new(
        cls: Type[T_TrackerOpeningProof],
        k_r_G: G1Point,
        r_G: G1Point,
        k_G: G1Point,
        k: Scalar,
        transcript: CurdleproofsTranscript,
    ) -> T_TrackerOpeningProof:
        blinder = generate_blinders(1)[0]
        A = G1 * blinder
        B = r_G * blinder

        transcript.append_list(
            b"tracker_opening_proof",
            points_projective_to_bytes([k_G, G1, k_r_G, r_G, A, B]),
        )

        challenge = transcript.get_and_append_challenge(
            b"tracker_opening_proof_challenge"
        )
        s = blinder - challenge * k

        return cls(A, B, s)

    def verify(
        self,
        transcript: CurdleproofsTranscript,
        k_r_G: G1Point,
        r_G: G1Point,
        k_G: G1Point,
    ):
        transcript.append_list(
            b"tracker_opening_proof",
            points_projective_to_bytes([k_G, G1, k_r_G, r_G, self.A, self.B]),
        )
        challenge = transcript.get_and_append_challenge(
            b"tracker_opening_proof_challenge"
        )

        Aprime = G1 * self.s + k_G * challenge
        Bprime = r_G * self.s + k_r_G * challenge

        assert Aprime == self.A and Bprime == self.B

    def to_json(self):
        return {
            "A": point_projective_to_json(self.A),
            "B": point_projective_to_json(self.B),
            "s": field_to_json(self.s),
        }

    @classmethod
    def from_json(cls: Type[T_TrackerOpeningProof], json) -> T_TrackerOpeningProof:
        return cls(
            point_projective_from_json(json["A"]),
            point_projective_from_json(json["B"]),
            field_from_json(json["s"]),
        )

    def to_bytes(self) -> bytes:
        return b''.join([
            g1_to_bytes(self.A),
            g1_to_bytes(self.B),
            fr_to_bytes(self.s)
        ])

    @classmethod
    def from_bytes(cls: Type[T_TrackerOpeningProof], b: BufReader) -> T_TrackerOpeningProof:
        return cls(
            A=b.read_g1(),
            B=b.read_g1(),
            s=b.read_fr()
        )
