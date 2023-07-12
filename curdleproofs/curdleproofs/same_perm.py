from functools import reduce
from curdleproofs.grand_prod import GrandProductProof
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Type, TypeVar
from curdleproofs.util import (
    field_to_bytes,
    point_projective_from_json,
    point_projective_to_json,
    points_projective_to_bytes,
    get_permutation,
    BufReader,
    g1_to_bytes,
)
from curdleproofs.msm_accumulator import MSMAccumulator, compute_MSM
from operator import mul as op_mul
from py_arkworks_bls12381 import G1Point, Scalar

T_SAME_PERM_PROOF = TypeVar("T_SAME_PERM_PROOF", bound="SamePermutationProof")


class SamePermutationProof:
    def __init__(self, B: G1Point, grand_prod_proof: GrandProductProof) -> None:
        self.B = B
        self.grand_prod_proof = grand_prod_proof

    @classmethod
    def new(
        cls: Type[T_SAME_PERM_PROOF],
        crs_G_vec: List[G1Point],
        crs_H_vec: List[G1Point],
        crs_U: G1Point,
        A: G1Point,
        M: G1Point,
        vec_a: List[Scalar],
        permutation: List[int],
        vec_a_blinders: List[Scalar],
        vec_m_blinders: List[Scalar],
        transcript: CurdleproofsTranscript,
    ) -> T_SAME_PERM_PROOF:
        n_blinders = len(vec_a_blinders)
        ell = len(crs_G_vec)

        transcript.append_list(b"same_perm_step1", points_projective_to_bytes([A, M]))
        transcript.append_list(b"same_perm_step1", list(map(field_to_bytes, vec_a)))
        alpha = transcript.get_and_append_challenge(b"same_perm_alpha")
        beta = transcript.get_and_append_challenge(b"same_perm_beta")

        vec_a_permuted = get_permutation(vec_a, permutation)
        permuted_polynomial_factors = [
            a + Scalar(m) * alpha + beta for (a, m) in zip(vec_a_permuted, permutation)
        ]
        gprod_result = reduce(op_mul, permuted_polynomial_factors, Scalar(1))

        vec_beta_repeated = [beta] * ell
        B = (A + M * alpha) + compute_MSM(crs_G_vec, vec_beta_repeated)

        vec_b_blinders = [
            vec_a_blinders[i] + alpha * vec_m_blinders[i] for i in range(0, n_blinders)
        ]

        grand_product_proof = GrandProductProof.new(
            crs_G_vec=crs_G_vec,
            crs_H_vec=crs_H_vec,
            crs_U=crs_U,
            B=B,
            gprod_result=gprod_result,
            vec_b=permuted_polynomial_factors,
            vec_b_blinders=vec_b_blinders,
            transcript=transcript,
        )

        return cls(B, grand_product_proof)

    def verify(
        self,
        crs_G_vec: List[G1Point],
        crs_H_vec: List[G1Point],
        crs_U: G1Point,
        crs_G_sum: G1Point,
        crs_H_sum: G1Point,
        A: G1Point,
        M: G1Point,
        vec_a: List[Scalar],
        n_blinders: int,
        transcript: CurdleproofsTranscript,
        msm_accumulator: MSMAccumulator,
    ):
        ell = len(crs_G_vec)

        # Step 1
        transcript.append_list(b"same_perm_step1", points_projective_to_bytes([A, M]))
        transcript.append_list(b"same_perm_step1", list(map(field_to_bytes, vec_a)))

        alpha = transcript.get_and_append_challenge(b"same_perm_alpha")
        beta = transcript.get_and_append_challenge(b"same_perm_beta")

        # Step 2
        polynomial_factors = [
            a + Scalar(i) * alpha + beta for (a, i) in zip(vec_a, range(0, ell))
        ]
        gprod_result = reduce(op_mul, polynomial_factors, Scalar(1))
        vec_beta_repeated = [beta] * ell
        msm_accumulator.accumulate_check(
            (self.B - A) - (M * alpha),
            crs_G_vec,
            vec_beta_repeated,
        )

        self.grand_prod_proof.verify(
            crs_G_vec=crs_G_vec,
            crs_H_vec=crs_H_vec,
            crs_U=crs_U,
            crs_G_sum=crs_G_sum,
            crs_H_sum=crs_H_sum,
            B=self.B,
            gprod_result=gprod_result,
            n_blinders=n_blinders,
            transcript=transcript,
            msm_accumulator=msm_accumulator,
        )

    def to_json(self):
        return {
            "B": point_projective_to_json(self.B),
            "grand_prod_proof": self.grand_prod_proof.to_json(),
        }

    @classmethod
    def from_json(cls: Type[T_SAME_PERM_PROOF], json) -> T_SAME_PERM_PROOF:
        return cls(
            B=point_projective_from_json(json["B"]),
            grand_prod_proof=GrandProductProof.from_json(json["grand_prod_proof"]),
        )

    def to_bytes(self) -> bytes:
        return b''.join([
            g1_to_bytes(self.B),
            self.grand_prod_proof.to_bytes(),
        ])

    @classmethod
    def from_bytes(cls: Type[T_SAME_PERM_PROOF], b: BufReader, n: int) -> T_SAME_PERM_PROOF:
        return cls(
            B=b.read_g1(),
            grand_prod_proof=GrandProductProof.from_bytes(b, n),
        )
