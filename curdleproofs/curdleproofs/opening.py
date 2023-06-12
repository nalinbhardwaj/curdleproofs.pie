import json
from typing import Type, TypeVar
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript

from curdleproofs.util import (
    Fr,
    PointProjective,
    field_from_json,
    field_to_json,
    generate_blinders,
    point_projective_from_json,
    point_projective_to_bytes,
    point_projective_to_json,
    points_projective_to_bytes,
)
from py_ecc.optimized_bls12_381.optimized_curve import (
    curve_order,
    G1,
    multiply,
    normalize,
    add,
    neg,
    eq,
)
from py_ecc.bls.g2_primitives import (
    G1_to_pubkey,
    pubkey_to_G1,
)
from py_ecc.bls.hash import i2osp, os2ip

T_TrackerOpeningProof = TypeVar("T_TrackerOpeningProof", bound="TrackerOpeningProof")


# proof of knowledge of `k` s.t. k_r_G = k*r_G and k_G = k*G
class TrackerOpeningProof:
    def __init__(
        self,
        A: PointProjective,
        B: PointProjective,
        s: Fr,
    ) -> None:
        self.A = A
        self.B = B
        self.s = s

    @classmethod
    def new(
        cls: Type[T_TrackerOpeningProof],
        k_r_G: PointProjective,
        r_G: PointProjective,
        k_G: PointProjective,
        G: PointProjective,
        k: Fr,
        transcript: CurdleproofsTranscript,
    ) -> T_TrackerOpeningProof:
        blinder = generate_blinders(1)[0]
        A = multiply(G, int(blinder))
        B = multiply(r_G, int(blinder))

        transcript.append_list(
            b"tracker_opening_proof",
            points_projective_to_bytes([k_G, G, k_r_G, r_G, A, B]),
        )

        challenge = transcript.get_and_append_challenge(
            b"tracker_opening_proof_challenge"
        )
        s = blinder - challenge * k

        return cls(A, B, s)

    def verify(
        self,
        transcript: CurdleproofsTranscript,
        k_r_G: PointProjective,
        r_G: PointProjective,
        k_G: PointProjective,
    ) -> bool:
        transcript.append_list(
            b"tracker_opening_proof",
            points_projective_to_bytes([k_G, G1, k_r_G, r_G, self.A, self.B]),
        )
        challenge = transcript.get_and_append_challenge(
            b"tracker_opening_proof_challenge"
        )

        Aprime = add(multiply(G1, int(self.s)), multiply(k_G, int(challenge)))
        Bprime = add(multiply(r_G, int(self.s)), multiply(k_r_G, int(challenge)))

        return eq(Aprime, self.A) and eq(Bprime, self.B)

    def to_bytes(self) -> bytes:
        A = G1_to_pubkey(self.A)
        B = G1_to_pubkey(self.B)
        s = i2osp(int(self.s), 32)
        return A + B + s
    
    @classmethod
    def from_bytes(cls: Type[T_TrackerOpeningProof], b: bytes) -> T_TrackerOpeningProof:
        A = pubkey_to_G1(b[0:48])
        B = pubkey_to_G1(b[48:96])
        s = Fr(os2ip(b[96:128]))
        return cls(A, B, s)

