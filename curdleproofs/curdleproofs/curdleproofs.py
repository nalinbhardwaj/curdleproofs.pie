import json
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.ipa import generate_blinders
from curdleproofs.util import (
    point_projective_from_json,
    point_projective_to_bytes,
    point_projective_to_json,
    points_projective_to_bytes,
    BufReader,
    g1_to_bytes,
    random_scalar,
    Z1,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Tuple, Type, TypeVar
from curdleproofs.util import get_permutation, g1_is_inf
from curdleproofs.msm_accumulator import MSMAccumulator, compute_MSM
from curdleproofs.same_perm import SamePermutationProof
from curdleproofs.same_msm import SameMSMProof
from curdleproofs.same_scalar import SameScalarProof
from curdleproofs.commitment import GroupCommitment
from py_arkworks_bls12381 import G1Point, Scalar

N_BLINDERS = 4

T_CurdleProofsProof = TypeVar("T_CurdleProofsProof", bound="CurdleProofsProof")


class CurdleProofsProof:
    def __init__(
        self,
        A: G1Point,
        cm_T: GroupCommitment,
        cm_U: GroupCommitment,
        R: G1Point,
        S: G1Point,
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
        vec_R: List[G1Point],
        vec_S: List[G1Point],
        vec_T: List[G1Point],
        vec_U: List[G1Point],
        M: G1Point,
        permutation: List[int],
        k: Scalar,
        vec_m_blinders: List[Scalar],
    ) -> T_CurdleProofsProof:
        ell = len(vec_R)

        transcript = CurdleproofsTranscript(b"curdleproofs")

        transcript.append_list(
            b"curdleproofs_step1", points_projective_to_bytes(vec_R + vec_S + vec_T + vec_U)
        )
        transcript.append(b"curdleproofs_step1", point_projective_to_bytes(M))
        vec_a = transcript.get_and_append_challenges(b"curdleproofs_vec_a", ell)

        vec_a_blinders = generate_blinders(N_BLINDERS - 2)
        vec_r_a_prime = vec_a_blinders + [Scalar(0), Scalar(0)]
        vec_a_permuted = get_permutation(vec_a, permutation)

        A = compute_MSM(crs.vec_G, vec_a_permuted) + compute_MSM(crs.vec_H, vec_r_a_prime)

        same_perm_proof = SamePermutationProof.new(
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

        r_t = random_scalar()
        r_u = random_scalar()
        R = compute_MSM(vec_R, vec_a)
        S = compute_MSM(vec_S, vec_a)

        cm_T: GroupCommitment = GroupCommitment.new(
            crs.G_t, crs.H, R * k, r_t
        )
        cm_U: GroupCommitment = GroupCommitment.new(
            crs.G_u, crs.H, S * k, r_u
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

        A_prime = A + cm_T.T_1 + cm_U.T_1

        vec_G_with_blinders = (
            crs.vec_G + crs.vec_H[: (N_BLINDERS - 2)] + [crs.G_t, crs.G_u]
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
        vec_R: List[G1Point],
        vec_S: List[G1Point],
        vec_T: List[G1Point],
        vec_U: List[G1Point],
        M: G1Point,
    ):
        ell = len(vec_R)

        transcript = CurdleproofsTranscript(b"curdleproofs")
        msm_accumulator = MSMAccumulator()

        if g1_is_inf(vec_T[0]):
            raise Exception("vec_T[0] is infinity")

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

        A_prime = self.A + self.cm_T.T_1 + self.cm_U.T_1

        vec_G_with_blinders = (
            crs.vec_G + crs.vec_H[: (N_BLINDERS - 2)] + [crs.G_t, crs.G_u]
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

        msm_accumulator.verify()

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
    vec_R: List[G1Point],
    vec_S: List[G1Point],
    permutation: List[int],
    k: Scalar,
) -> Tuple[List[G1Point], List[G1Point], G1Point, List[Scalar]]:
    ell = len(crs.vec_G)

    vec_T = [R * k for R in vec_R]
    vec_U = [S * k for S in vec_S]
    vec_T = get_permutation(vec_T, permutation)
    vec_U = get_permutation(vec_U, permutation)

    range_as_fr = [Scalar(i) for i in range(ell)]
    sigma_ell = get_permutation(range_as_fr, permutation)

    vec_m_blinders = generate_blinders(N_BLINDERS)
    M = compute_MSM(crs.vec_G, sigma_ell) + compute_MSM(crs.vec_H, vec_m_blinders)

    return vec_T, vec_U, M, vec_m_blinders


T_VerifierInput = TypeVar("T_VerifierInput", bound="VerifierInput")


class VerifierInput:
    def __init__(
        self,
        vec_R: List[G1Point],
        vec_S: List[G1Point],
        vec_T: List[G1Point],
        vec_U: List[G1Point],
        M: G1Point,
    ) -> None:
        self.vec_R = vec_R
        self.vec_S = vec_S
        self.vec_T = vec_T
        self.vec_U = vec_U
        self.M = M

    def to_json(self) -> str:
        dic = {
            "vec_R": [point_projective_to_json(R) for R in self.vec_R],
            "vec_S": [point_projective_to_json(S) for S in self.vec_S],
            "vec_T": [point_projective_to_json(T) for T in self.vec_T],
            "vec_U": [point_projective_to_json(U) for U in self.vec_U],
            "M": point_projective_to_json(self.M),
        }
        return json.dumps(dic)

    @classmethod
    def from_json(cls: Type[T_VerifierInput], json_str: str) -> T_VerifierInput:
        dic = json.loads(json_str)
        return cls(
            vec_R=[point_projective_from_json(R) for R in dic["vec_R"]],
            vec_S=[point_projective_from_json(S) for S in dic["vec_S"]],
            vec_T=[point_projective_from_json(T) for T in dic["vec_T"]],
            vec_U=[point_projective_from_json(U) for U in dic["vec_U"]],
            M=point_projective_from_json(dic["M"]),
        )
