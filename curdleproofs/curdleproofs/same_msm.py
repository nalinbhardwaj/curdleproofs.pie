import json
from curdleproofs.util import (
    field_from_json,
    point_projective_from_json,
    points_projective_to_bytes,
    get_verification_scalars_bitstring,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Optional, Tuple, Type, TypeVar
from curdleproofs.util import PointProjective, Fr, invert
from curdleproofs.msm_accumulator import MSMAccumulator, compute_MSM
from py_ecc.optimized_bls12_381.optimized_curve import (
    multiply,
    add,
)

T_SameMSMProof = TypeVar("T_SameMSMProof", bound="SameMSMProof")


class SameMSMProof:
    def __init__(
        self,
        B_a: PointProjective,
        B_t: PointProjective,
        B_u: PointProjective,
        vec_L_A: List[PointProjective],
        vec_L_T: List[PointProjective],
        vec_L_U: List[PointProjective],
        vec_R_A: List[PointProjective],
        vec_R_T: List[PointProjective],
        vec_R_U: List[PointProjective],
        x_final: Fr,
    ) -> None:
        self.B_a = B_a
        self.B_t = B_t
        self.B_u = B_u
        self.vec_L_A = vec_L_A
        self.vec_L_T = vec_L_T
        self.vec_L_U = vec_L_U
        self.vec_R_A = vec_R_A
        self.vec_R_T = vec_R_T
        self.vec_R_U = vec_R_U
        self.x_final = x_final

    def verification_scalars(
        self, n: int, transcript: CurdleproofsTranscript
    ) -> Tuple[Optional[Tuple[List[Fr], List[Fr], List[Fr]]], str]:
        lg_n = len(self.vec_L_A)
        if lg_n >= 32:
            return None, "lg_n >= 32"
        if 2**lg_n != n:
            return None, "2**lg_n != n"

        bitstring = get_verification_scalars_bitstring(n, lg_n)

        challenges: List[Fr] = []
        for i in range(0, lg_n):
            transcript.append_list(
                b"same_msm_loop",
                points_projective_to_bytes(
                    [
                        self.vec_L_A[i],
                        self.vec_L_T[i],
                        self.vec_L_U[i],
                        self.vec_R_A[i],
                        self.vec_R_T[i],
                        self.vec_R_U[i],
                    ]
                ),
            )
            challenges.append(transcript.get_and_append_challenge(b"same_msm_gamma"))

        challenges_inv = list(map(invert, challenges))

        vec_s: List[Fr] = []
        for i in range(0, n):
            vec_s.append(Fr.one())
            for j in bitstring[i]:
                vec_s[i] *= challenges[j]

        return (challenges, challenges_inv, vec_s), ""

    def verify(
        self,
        crs_G_vec: List[PointProjective],
        A: PointProjective,
        Z_t: PointProjective,
        Z_u: PointProjective,
        vec_T: List[PointProjective],
        vec_U: List[PointProjective],
        transcript: CurdleproofsTranscript,
        msm_accumulator: MSMAccumulator,
    ) -> Tuple[bool, str]:
        n = len(vec_T)

        transcript.append_list(
            b"same_msm_step1", points_projective_to_bytes([A, Z_t, Z_u])
        )
        transcript.append_list(
            b"same_msm_step1", points_projective_to_bytes(vec_T + vec_U)
        )
        transcript.append_list(
            b"same_msm_step1",
            points_projective_to_bytes([self.B_a, self.B_t, self.B_u]),
        )
        alpha = transcript.get_and_append_challenge(b"same_msm_alpha")

        (ret, err) = self.verification_scalars(n, transcript)

        if ret is None:
            return False, err

        vec_gamma, vec_gamma_inv, vec_s = ret

        vec_x_times_s = [self.x_final * s_i for s_i in vec_s]

        A_a = add(self.B_a, multiply(A, int(alpha)))
        Z_t_a = add(self.B_t, multiply(Z_t, int(alpha)))
        Z_u_a = add(self.B_u, multiply(Z_u, int(alpha)))

        point_lhs = add(
            add(compute_MSM(self.vec_L_A, vec_gamma), A_a),
            compute_MSM(self.vec_R_A, vec_gamma_inv),
        )
        msm_accumulator.accumulate_check(point_lhs, crs_G_vec, vec_x_times_s)

        point_lhs = add(
            add(compute_MSM(self.vec_L_T, vec_gamma), Z_t_a),
            compute_MSM(self.vec_R_T, vec_gamma_inv),
        )
        msm_accumulator.accumulate_check(point_lhs, vec_T, vec_x_times_s)

        point_lhs = add(
            add(compute_MSM(self.vec_L_U, vec_gamma), Z_u_a),
            compute_MSM(self.vec_R_U, vec_gamma_inv),
        )
        msm_accumulator.accumulate_check(point_lhs, vec_U, vec_x_times_s)

        return True, ""

    @classmethod
    def from_json(cls: Type[T_SameMSMProof], json_str: str) -> T_SameMSMProof:
        dic = json.loads(json_str)
        return cls(
            B_a=point_projective_from_json(dic["B_a"]),
            B_t=point_projective_from_json(dic["B_t"]),
            B_u=point_projective_from_json(dic["B_u"]),
            vec_L_A=[point_projective_from_json(L_A) for L_A in dic["vec_L_A"]],
            vec_L_T=[point_projective_from_json(L_T) for L_T in dic["vec_L_T"]],
            vec_L_U=[point_projective_from_json(L_U) for L_U in dic["vec_L_U"]],
            vec_R_A=[point_projective_from_json(R_A) for R_A in dic["vec_R_A"]],
            vec_R_T=[point_projective_from_json(R_T) for R_T in dic["vec_R_T"]],
            vec_R_U=[point_projective_from_json(R_U) for R_U in dic["vec_R_U"]],
            x_final=field_from_json(dic["x_final"], Fr),
        )
