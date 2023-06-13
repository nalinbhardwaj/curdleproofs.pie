import json
from math import log2
import random
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.ipa import generate_blinders
from curdleproofs.util import (
    affine_to_projective,
    point_affine_to_bytes,
    point_projective_from_json,
    point_projective_to_bytes,
    point_projective_to_json,
    points_projective_to_bytes,
    get_random_point,
    BufReader,
    g1_to_bytes,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Optional, Tuple, Type, TypeVar
from curdleproofs.util import (
    PointProjective,
    Fr,
    field_to_bytes,
    invert,
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
    is_inf,
    FQ,
)
from curdleproofs.same_perm import SamePermutationProof
from curdleproofs.same_msm import SameMSMProof
from curdleproofs.same_scalar import SameScalarProof
from curdleproofs.commitment import GroupCommitment
from py_ecc.bls.g2_primitives import G1_to_pubkey, pubkey_to_G1
from eth_typing import BLSPubkey

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

    @classmethod
    def new(
        cls: Type[T_CurdleProofsProof],
        crs: CurdleproofsCrs,
        vec_R: List[PointProjective],
        vec_S: List[PointProjective],
        vec_T: List[PointProjective],
        vec_U: List[PointProjective],
        M: PointProjective,
        permutation: List[int],
        k: Fr,
        vec_m_blinders: List[Fr],
    ) -> T_CurdleProofsProof:
        ell = len(vec_R)

        transcript = CurdleproofsTranscript(b"curdleproofs")

        transcript.append_list(
            b"curdleproofs_step1", points_projective_to_bytes(vec_R + vec_S + vec_T + vec_U)
        )
        transcript.append(b"curdleproofs_step1", point_projective_to_bytes(M))
        vec_a = transcript.get_and_append_challenges(b"curdleproofs_vec_a", ell)

        vec_a_blinders = generate_blinders(crs.n_blinders() - 2)
        vec_r_a_prime = vec_a_blinders + [Fr.zero(), Fr.zero()]
        vec_a_permuted = get_permutation(vec_a, permutation)

        A = add(
            compute_MSM(crs.vec_G, vec_a_permuted),
            compute_MSM(crs.vec_H, vec_r_a_prime),
        )

        (same_perm_proof, err) = SamePermutationProof.new(
            crs_G_vec=crs.vec_G,
            crs_H_vec=crs.vec_H,
            crs_U=crs.H,
            A=A,
            M=M,
            vec_a=vec_a,
            permutation=permutation,
            vec_a_blinders=vec_r_a_prime,
            vec_m_blinders=vec_m_blinders,
            transcript=transcript,
        )

        if same_perm_proof is None:
            raise Exception(err)

        r_t = Fr(random.randint(1, Fr.field_modulus))
        r_u = Fr(random.randint(1, Fr.field_modulus))
        R = compute_MSM(vec_R, vec_a)
        S = compute_MSM(vec_S, vec_a)

        cm_T: GroupCommitment = GroupCommitment.new(
            crs.G_t, crs.H, multiply(R, int(k)), r_t
        )
        cm_U: GroupCommitment = GroupCommitment.new(
            crs.G_u, crs.H, multiply(S, int(k)), r_u
        )

        same_scalar_proof = SameScalarProof.new(
            crs_G_t=crs.G_t,
            crs_G_u=crs.G_u,
            crs_H=crs.H,
            R=R,
            S=S,
            cm_T=cm_T,
            cm_U=cm_U,
            k=k,
            r_t=r_t,
            r_u=r_u,
            transcript=transcript,
        )

        A_prime = add(add(A, cm_T.T_1), cm_U.T_1)

        vec_G_with_blinders = (
            crs.vec_G + crs.vec_H[: (crs.n_blinders() - 2)] + [crs.G_t, crs.G_u]
        )

        vec_T_with_blinders = vec_T + [
            Z1,
            Z1,
            crs.H,
            Z1,
        ]

        vec_U_with_blinders = vec_U + [
            Z1,
            Z1,
            Z1,
            crs.H,
        ]

        vec_a_with_blinders = vec_a_permuted + vec_a_blinders + [r_t, r_u]

        same_msm_proof = SameMSMProof.new(
            crs_G_vec=vec_G_with_blinders,
            A=A_prime,
            Z_t=cm_T.T_2,
            Z_u=cm_U.T_2,
            vec_T=vec_T_with_blinders,
            vec_U=vec_U_with_blinders,
            vec_x=vec_a_with_blinders,
            transcript=transcript,
        )

        return cls(
            A=A,
            cm_T=cm_T,
            cm_U=cm_U,
            R=R,
            S=S,
            same_perm_proof=same_perm_proof,
            same_scalar_proof=same_scalar_proof,
            same_msm_proof=same_msm_proof,
        )

    def verify(
        self,
        crs: CurdleproofsCrs,
        vec_R: List[PointProjective],
        vec_S: List[PointProjective],
        vec_T: List[PointProjective],
        vec_U: List[PointProjective],
        M: PointProjective,
    ) -> Tuple[bool, str]:
        ell = len(vec_R)

        transcript = CurdleproofsTranscript(b"curdleproofs")
        msm_accumulator = MSMAccumulator()

        if is_inf(vec_T[0]):
            return False, "vec_T[0] is infinity"

        transcript.append_list(
            b"curdleproofs_step1", points_projective_to_bytes(vec_R + vec_S + vec_T + vec_U)
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
            n_blinders=crs.n_blinders(),
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
            crs.vec_G + crs.vec_H[: (crs.n_blinders() - 2)] + [crs.G_t, crs.G_u]
        )

        vec_T_with_blinders = vec_T + [
            Z1,
            Z1,
            crs.H,
            Z1,
        ]

        vec_U_with_blinders = vec_U + [
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
            self.R, vec_R, vec_a
        )
        msm_accumulator.accumulate_check(
            self.S, vec_S, vec_a
        )

        msm_verify = msm_accumulator.verify()

        if not msm_verify:
            return False, "MSM check failed"

        return True, ""

    def to_json(self):
        return {
            "A": point_projective_to_json(self.A),
            "cm_T": self.cm_T.to_json(),
            "cm_U": self.cm_U.to_json(),
            "R": point_projective_to_json(self.R),
            "S": point_projective_to_json(self.S),
            "same_perm_proof": self.same_perm_proof.to_json(),
            "same_scalar_proof": self.same_scalar_proof.to_json(),
            "same_msm_proof": self.same_msm_proof.to_json(),
        }

    @classmethod
    def from_json(cls: Type[T_CurdleProofsProof], json) -> T_CurdleProofsProof:
        return cls(
            A=point_projective_from_json(json["A"]),
            cm_T=GroupCommitment.from_json(json["cm_T"]),
            cm_U=GroupCommitment.from_json(json["cm_U"]),
            R=point_projective_from_json(json["R"]),
            S=point_projective_from_json(json["S"]),
            same_perm_proof=SamePermutationProof.from_json(json["same_perm_proof"]),
            same_scalar_proof=SameScalarProof.from_json(json["same_scalar_proof"]),
            same_msm_proof=SameMSMProof.from_json(json["same_msm_proof"]),
        )
    
    def to_bytes(self) -> bytes:
        return b''.join([
            g1_to_bytes(self.A),
            self.cm_T.to_bytes(),
            self.cm_U.to_bytes(),
            g1_to_bytes(self.R),
            g1_to_bytes(self.S),
            self.same_perm_proof.to_bytes(),
            self.same_scalar_proof.to_bytes(),
            self.same_msm_proof.to_bytes(),
        ])

    @classmethod
    def from_bytes(cls: Type[T_CurdleProofsProof], b: BufReader, n: int) -> T_CurdleProofsProof:
        return cls(
            A=b.read_g1(),
            cm_T=GroupCommitment.from_bytes(b),
            cm_U=GroupCommitment.from_bytes(b),
            R=b.read_g1(),
            S=b.read_g1(),
            same_perm_proof=SamePermutationProof.from_bytes(b, n),
            same_scalar_proof=SameScalarProof.from_bytes(b),
            same_msm_proof=SameMSMProof.from_bytes(b, n),
        )


def shuffle_permute_and_commit_input(
    crs: CurdleproofsCrs,
    vec_R: List[PointProjective],
    vec_S: List[PointProjective],
    permutation: List[int],
    k: Fr,
) -> Tuple[List[PointProjective], List[PointProjective], PointProjective, List[Fr]]:
    ell = len(crs.vec_G)

    vec_T = [multiply(R, int(k)) for R in vec_R]
    vec_U = [multiply(S, int(k)) for S in vec_S]
    vec_T = get_permutation(vec_T, permutation)
    vec_U = get_permutation(vec_U, permutation)

    range_as_fr = [Fr(i) for i in range(ell)]
    sigma_ell = get_permutation(range_as_fr, permutation)

    vec_m_blinders = generate_blinders(crs.n_blinders())
    M = add(compute_MSM(crs.vec_G, sigma_ell), compute_MSM(crs.vec_H, vec_m_blinders))

    return vec_T, vec_U, M, vec_m_blinders
