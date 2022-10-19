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

T_TrackerOpeningProof = TypeVar("T_TrackerOpeningProof", bound="TrackerOpeningProof")

# proof of knowledge of `k` s.t. k_r_G = k*r_G and k_G = k*G
class TrackerOpeningProof:
    def __init__(
        self,
        k_r_G: PointProjective,
        r_G: PointProjective,
        k_G: PointProjective,
        G: PointProjective,
        A: PointProjective,
        B: PointProjective,
        s: Fr,
    ) -> None:
        self.k_r_G = k_r_G
        self.r_G = r_G
        self.k_G = k_G
        self.G = G
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

        return cls(k_r_G, r_G, k_G, G, A, B, s)

    def verify(
        self,
        transcript: CurdleproofsTranscript,
    ) -> bool:
        transcript.append_list(
            b"tracker_opening_proof",
            points_projective_to_bytes(
                [self.k_G, self.G, self.k_r_G, self.r_G, self.A, self.B]
            ),
        )
        challenge = transcript.get_and_append_challenge(
            b"tracker_opening_proof_challenge"
        )

        Aprime = add(multiply(self.G, int(self.s)), multiply(self.k_G, int(challenge)))
        Bprime = add(
            multiply(self.r_G, int(self.s)), multiply(self.k_r_G, int(challenge))
        )

        return eq(Aprime, self.A) and eq(Bprime, self.B)

    def to_json(self) -> str:
        return json.dumps(
            {
                "k_r_G": point_projective_to_json(self.k_r_G),
                "r_G": point_projective_to_json(self.r_G),
                "k_G": point_projective_to_json(self.k_G),
                "G": point_projective_to_json(self.G),
                "A": point_projective_to_json(self.A),
                "B": point_projective_to_json(self.B),
                "s": field_to_json(self.s),
            }
        )

    @classmethod
    def from_json(
        cls: Type[T_TrackerOpeningProof], json_str: str
    ) -> T_TrackerOpeningProof:
        dic = json.loads(json_str)
        return cls(
            point_projective_from_json(dic["k_r_G"]),
            point_projective_from_json(dic["r_G"]),
            point_projective_from_json(dic["k_G"]),
            point_projective_from_json(dic["G"]),
            point_projective_from_json(dic["A"]),
            point_projective_from_json(dic["B"]),
            field_from_json(dic["s"], Fr),
        )
