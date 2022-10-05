from functools import reduce
import json

from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.grand_prod import GrandProductProof
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Optional, Tuple, Type, TypeVar
from curdleproofs.util import (
    PointAffine,
    PointProjective,
    Fr,
    field_to_bytes,
    point_projective_from_json,
    point_projective_to_json,
    points_projective_to_bytes,
    get_random_point,
    get_permutation,
)
from curdleproofs.msm_accumulator import MSMAccumulator, compute_MSM
from py_ecc.optimized_bls12_381.optimized_curve import (
    curve_order,
    G1,
    multiply,
    normalize,
    add,
    neg,
    Z1,
)
from operator import mul as op_mul

T_SAME_PERM_PROOF = TypeVar("T_SAME_PERM_PROOF", bound="SamePermutationProof")


class SamePermutationProof:
    def __init__(self, B: PointProjective, grand_prod_proof: GrandProductProof) -> None:
        self.B = B
        self.grand_prod_proof = grand_prod_proof

    def verify(
        self,
        crs_G_vec: List[PointProjective],
        crs_H_vec: List[PointProjective],
        crs_U: PointProjective,
        crs_G_sum: PointProjective,
        crs_H_sum: PointProjective,
        A: PointProjective,
        M: PointProjective,
        vec_a: List[Fr],
        n_blinders: int,
        transcript: CurdleproofsTranscript,
        msm_accumulator: MSMAccumulator,
    ) -> Tuple[bool, str]:
        ell = len(crs_G_vec)

        # Step 1
        transcript.append_list(b"same_perm_step1", points_projective_to_bytes([A, M]))
        transcript.append_list(b"same_perm_step1", list(map(field_to_bytes, vec_a)))

        alpha = transcript.get_and_append_challenge(b"same_perm_alpha")
        beta = transcript.get_and_append_challenge(b"same_perm_beta")

        # Step 2
        polynomial_factors = [
            a + Fr(i) * alpha + beta for (a, i) in zip(vec_a, range(0, ell))
        ]
        gprod_result = reduce(op_mul, polynomial_factors, Fr.one())
        vec_beta_repeated = [beta] * ell
        msm_accumulator.accumulate_check(
            add(add(self.B, neg(A)), neg(multiply(M, int(alpha)))),
            crs_G_vec,
            vec_beta_repeated,
        )

        (grand_prod_verify, err) = self.grand_prod_proof.verify(
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

        if not grand_prod_verify:
            return (False, err)

        return (True, "")

    @classmethod
    def from_json(cls: Type[T_SAME_PERM_PROOF], json_str: str) -> T_SAME_PERM_PROOF:
        dic = json.loads(json_str)
        return cls(
            B=point_projective_from_json(dic["B"]),
            grand_prod_proof=GrandProductProof.from_json(dic["grand_prod_proof"]),
        )
