import json
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.util import (
    affine_to_projective,
    point_affine_from_json,
    point_projective_from_json,
    point_projective_to_bytes,
    points_affine_to_bytes,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Optional, Tuple, Type, TypeVar
from curdleproofs.util import PointAffine, PointProjective
from curdleproofs.msm_accumulator import MSMAccumulator, compute_MSM
from py_ecc.optimized_bls12_381.optimized_curve import (
    curve_order,
    G1,
    multiply,
    normalize,
    add,
    neg,
    Z1,
    is_inf,
    FQ,
)
from curdleproofs.same_perm import SamePermutationProof
from curdleproofs.same_msm import SameMSMProof
from curdleproofs.same_scalar import SameScalarProof
from curdleproofs.commitment import GroupCommitment

N_BLINDERS = 4

T_CurdleProofsProof = TypeVar("T_CurdleProofsProof", bound="CurdleProofsProof")


class CurdleProofsProof:
    def __init__(
        self,
        A: PointProjective,
        cm_T: GroupCommitment,
        cm_U: GroupCommitment,
        R: PointProjective,
        S: PointProjective,
        same_perm_proof: SamePermutationProof,
        same_scalar_proof: SameScalarProof,
        same_msm_proof: SameMSMProof,
    ) -> None:
        self.A = A
        self.cm_T = cm_T
        self.cm_U = cm_U
        self.R = R
        self.S = S
        self.same_perm_proof = same_perm_proof
        self.same_scalar_proof = same_scalar_proof
        self.same_msm_proof = same_msm_proof

    def verify(
        self,
        crs: CurdleproofsCrs,
        vec_R: List[PointAffine],
        vec_S: List[PointAffine],
        vec_T: List[PointAffine],
        vec_U: List[PointAffine],
        M: PointProjective,
    ) -> Tuple[bool, str]:
        ell = len(vec_R)

        transcript = CurdleproofsTranscript(b"curdleproofs")
        msm_accumulator = MSMAccumulator()

        if is_inf(affine_to_projective(vec_T[0])):
            return False, "vec_T[0] is infinity"

        transcript.append_list(
            b"curdleproofs_step1", points_affine_to_bytes(vec_R + vec_S + vec_T + vec_U)
        )
        transcript.append(b"curdleproofs_step1", point_projective_to_bytes(M))
        vec_a = transcript.get_and_append_challenges(b"curdleproofs_vec_a", ell)

        self.same_perm_proof.verify(
            crs_G_vec=crs.vec_G,
            crs_H_vec=crs.vec_H,
            crs_U=crs.H,
            crs_G_sum=crs.G_sum,
            crs_H_sum=crs.H_sum,
            A=self.A,
            M=M,
            vec_a=vec_a,
            n_blinders=N_BLINDERS,
            transcript=transcript,
            msm_accumulator=msm_accumulator,
        )

        self.same_scalar_proof.verify(
            crs_G_t=crs.G_t,
            crs_G_u=crs.G_u,
            crs_H=crs.H,
            R=self.R,
            S=self.S,
            cm_T=self.cm_T,
            cm_U=self.cm_U,
            transcript=transcript,
        )

        A_prime = add(add(self.A, self.cm_T.T_1), self.cm_U.T_1)

        vec_G_with_blinders = (
            crs.vec_G + crs.vec_H[: (N_BLINDERS - 2)] + [crs.G_t, crs.G_u]
        )

        vec_T_with_blinders = list(map(affine_to_projective, vec_T)) + [
            Z1,
            Z1,
            crs.H,
            Z1,
        ]

        vec_U_with_blinders = list(map(affine_to_projective, vec_U)) + [
            Z1,
            Z1,
            Z1,
            crs.H,
        ]

        self.same_msm_proof.verify(
            crs_G_vec=vec_G_with_blinders,
            A=A_prime,
            Z_t=self.cm_T.T_2,
            Z_u=self.cm_U.T_2,
            vec_T=vec_T_with_blinders,
            vec_U=vec_U_with_blinders,
            transcript=transcript,
            msm_accumulator=msm_accumulator,
        )

        msm_accumulator.accumulate_check(
            self.R, list(map(affine_to_projective, vec_R)), vec_a
        )
        msm_accumulator.accumulate_check(
            self.S, list(map(affine_to_projective, vec_S)), vec_a
        )

        msm_verify = msm_accumulator.verify()

        if not msm_verify:
            return False, "MSM check failed"

        return True, ""

    @classmethod
    def from_json(cls: Type[T_CurdleProofsProof], json_str: str) -> T_CurdleProofsProof:
        dic = json.loads(json_str)
        return cls(
            A=point_projective_from_json(dic["A"]),
            cm_T=GroupCommitment.from_json(dic["cm_T"]),
            cm_U=GroupCommitment.from_json(dic["cm_U"]),
            R=point_projective_from_json(dic["R"]),
            S=point_projective_from_json(dic["S"]),
            same_perm_proof=SamePermutationProof.from_json(dic["same_perm_proof"]),
            same_scalar_proof=SameScalarProof.from_json(dic["same_scalar_proof"]),
            same_msm_proof=SameMSMProof.from_json(dic["same_msm_proof"]),
        )


T_VerifierInput = TypeVar("T_VerifierInput", bound="VerifierInput")


class VerifierInput:
    def __init__(
        self,
        vec_R: List[PointAffine],
        vec_S: List[PointAffine],
        vec_T: List[PointAffine],
        vec_U: List[PointAffine],
        M: PointProjective,
    ) -> None:
        self.vec_R = vec_R
        self.vec_S = vec_S
        self.vec_T = vec_T
        self.vec_U = vec_U
        self.M = M

    @classmethod
    def from_json(cls: Type[T_VerifierInput], json_str: str) -> T_VerifierInput:
        dic = json.loads(json_str)
        return cls(
            vec_R=[point_affine_from_json(R) for R in dic["vec_R"]],
            vec_S=[point_affine_from_json(S) for S in dic["vec_S"]],
            vec_T=[point_affine_from_json(T) for T in dic["vec_T"]],
            vec_U=[point_affine_from_json(U) for U in dic["vec_U"]],
            M=point_projective_from_json(dic["M"]),
        )
