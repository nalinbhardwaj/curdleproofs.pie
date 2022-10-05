from functools import reduce
import json
import operator
import random
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.ipa import IPA
from curdleproofs.util import (
    field_from_json,
    field_to_json,
    invert,
    point_projective_from_json,
    point_projective_to_bytes,
    get_random_point,
    point_projective_to_json,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Optional, Tuple, TypeVar, Type
from curdleproofs.util import (
    PointAffine,
    PointProjective,
    Fr,
    field_to_bytes,
    affine_to_projective,
)
from curdleproofs.msm_accumulator import MSMAccumulator, compute_MSM
from py_ecc.optimized_bls12_381.optimized_curve import (
    curve_order,
    G1,
    multiply,
    normalize,
    add,
    neg,
    eq,
    Z1,
)

T_GrandProductProof = TypeVar("T_GrandProductProof", bound="GrandProductProof")


class GrandProductProof:
    def __init__(self, C: PointProjective, r_p: Fr, ipa_proof: IPA):
        self.C = C
        self.r_p = r_p
        self.ipa_proof = ipa_proof

    def verify(
        self,
        crs_G_vec: List[PointProjective],
        crs_H_vec: List[PointProjective],
        crs_U: PointProjective,
        crs_G_sum: PointProjective,
        crs_H_sum: PointProjective,
        B: PointProjective,
        gprod_result: Fr,
        n_blinders: int,
        transcript: CurdleproofsTranscript,
        msm_accumulator: MSMAccumulator,
    ) -> Tuple[bool, str]:
        ell = len(crs_G_vec)

        # Step 1
        transcript.append(b"gprod_step1", point_projective_to_bytes(B))
        transcript.append(b"gprod_step1", field_to_bytes(gprod_result))
        alpha = transcript.get_and_append_challenge(b"gprod_alpha")

        # Step 2
        transcript.append(b"gprod_step2", point_projective_to_bytes(self.C))
        transcript.append(b"gprod_step2", field_to_bytes(self.r_p))
        beta = transcript.get_and_append_challenge(b"gprod_beta")
        beta_inv = invert(beta)

        # Step 3
        # Build `vec_u` for the optimization trick
        vec_u: List[Fr] = []
        pow_beta_inv = beta_inv
        for _ in range(0, ell):
            vec_u.append(pow_beta_inv)
            pow_beta_inv *= beta_inv

        vec_u.extend([beta_inv ** (ell + 1) for _ in range(0, n_blinders)])

        # Compute D
        D = add(
            add(B, neg(multiply(crs_G_sum, int(beta_inv)))),
            multiply(crs_H_sum, int(alpha)),
        )

        # Step 4
        # Build G
        vec_G = crs_G_vec + crs_H_vec

        inner_prod = (
            self.r_p * (beta ** (ell + 1)) + gprod_result * (beta**ell) - Fr.one()
        )

        (ipa_result, err) = self.ipa_proof.verify(
            crs_G_vec=vec_G,
            crs_H=crs_U,
            C=self.C,
            D=D,
            inner_prod=inner_prod,
            vec_u=vec_u,
            transcript=transcript,
            msm_accumulator=msm_accumulator,
        )

        if not ipa_result:
            return False, err

        return True, ""

    @classmethod
    def from_json(cls: Type[T_GrandProductProof], json_str: str) -> T_GrandProductProof:
        dic = json.loads(json_str)
        return cls(
            C=point_projective_from_json(dic["C"]),
            r_p=field_from_json(dic["r_p"], Fr),
            ipa_proof=IPA.from_json(dic["ipa_proof"]),
        )
