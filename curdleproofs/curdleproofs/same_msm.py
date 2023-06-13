import json
from math import log2
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.util import (
    affine_to_projective,
    field_from_json,
    field_to_json,
    point_affine_to_bytes,
    point_projective_from_json,
    point_projective_to_bytes,
    point_projective_to_json,
    points_affine_to_bytes,
    points_projective_to_bytes,
    get_random_point,
    generate_blinders,
    get_verification_scalars_bitstring,
    g1_from_bytes,
    BufReader,
    g1_to_bytes,
    fr_to_bytes,
    g1_list_to_bytes,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Optional, Tuple, Type, TypeVar
from curdleproofs.util import PointProjective, Fr, field_to_bytes, invert
from curdleproofs.msm_accumulator import MSMAccumulator, compute_MSM
from py_ecc.optimized_bls12_381.optimized_curve import (
    curve_order,
    G1,
    multiply,
    normalize,
    add,
    neg,
)
from py_ecc.bls.g2_primitives import G1_to_pubkey, pubkey_to_G1
from eth_typing import BLSPubkey

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

    @classmethod
    def new(
        cls: Type[T_SameMSMProof],
        crs_G_vec: List[PointProjective],
        A: PointProjective,
        Z_t: PointProjective,
        Z_u: PointProjective,
        vec_T: List[PointProjective],
        vec_U: List[PointProjective],
        vec_x: List[Fr],
        transcript: CurdleproofsTranscript,
    ) -> T_SameMSMProof:
        n = len(vec_x)
        lg_n = int(log2(n))
        assert 2**lg_n == n

        vec_L_T: List[PointProjective] = []
        vec_R_T: List[PointProjective] = []
        vec_L_U: List[PointProjective] = []
        vec_R_U: List[PointProjective] = []
        vec_L_A: List[PointProjective] = []
        vec_R_A: List[PointProjective] = []

        vec_r = generate_blinders(n)

        B_a = compute_MSM(crs_G_vec, vec_r)
        B_t = compute_MSM(vec_T, vec_r)
        B_u = compute_MSM(vec_U, vec_r)

        transcript.append_list(
            b"same_msm_step1", points_projective_to_bytes([A, Z_t, Z_u])
        )
        transcript.append_list(
            b"same_msm_step1", points_projective_to_bytes(vec_T + vec_U)
        )
        transcript.append_list(
            b"same_msm_step1", points_projective_to_bytes([B_a, B_t, B_u])
        )
        alpha = transcript.get_and_append_challenge(b"same_msm_alpha")

        for i in range(0, n):
            vec_x[i] = vec_r[i] + (alpha * vec_x[i])

        while len(vec_x) > 1:
            n //= 2

            x_L, x_R = vec_x[:n], vec_x[n:]
            T_L, T_R = vec_T[:n], vec_T[n:]
            U_L, U_R = vec_U[:n], vec_U[n:]
            G_L, G_R = crs_G_vec[:n], crs_G_vec[n:]

            L_A = compute_MSM(G_R, x_L)
            L_T = compute_MSM(T_R, x_L)
            L_U = compute_MSM(U_R, x_L)
            R_A = compute_MSM(G_L, x_R)
            R_T = compute_MSM(T_L, x_R)
            R_U = compute_MSM(U_L, x_R)

            vec_L_A.append(L_A)
            vec_L_T.append(L_T)
            vec_L_U.append(L_U)
            vec_R_A.append(R_A)
            vec_R_T.append(R_T)
            vec_R_U.append(R_U)

            transcript.append_list(
                b"same_msm_loop",
                points_projective_to_bytes([L_A, L_T, L_U, R_A, R_T, R_U]),
            )
            gamma = transcript.get_and_append_challenge(b"same_msm_gamma")
            gamma_inv = invert(gamma)

            for i in range(0, n):
                x_L[i] += gamma_inv * x_R[i]
                T_L[i] = add(T_L[i], multiply(T_R[i], int(gamma)))
                U_L[i] = add(U_L[i], multiply(U_R[i], int(gamma)))
                G_L[i] = add(G_L[i], multiply(G_R[i], int(gamma)))

            vec_x = x_L
            vec_T = T_L
            vec_U = U_L
            crs_G_vec = G_L

        return cls(
            B_a=B_a,
            B_t=B_t,
            B_u=B_u,
            vec_L_A=vec_L_A,
            vec_L_T=vec_L_T,
            vec_L_U=vec_L_U,
            vec_R_A=vec_R_A,
            vec_R_T=vec_R_T,
            vec_R_U=vec_R_U,
            x_final=vec_x[0],
        )

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

    def to_json(self) -> str:
        dic = {
            "B_a": point_projective_to_json(self.B_a),
            "B_t": point_projective_to_json(self.B_t),
            "B_u": point_projective_to_json(self.B_u),
            "vec_L_A": [point_projective_to_json(L_A) for L_A in self.vec_L_A],
            "vec_L_T": [point_projective_to_json(L_T) for L_T in self.vec_L_T],
            "vec_L_U": [point_projective_to_json(L_U) for L_U in self.vec_L_U],
            "vec_R_A": [point_projective_to_json(R_A) for R_A in self.vec_R_A],
            "vec_R_T": [point_projective_to_json(R_T) for R_T in self.vec_R_T],
            "vec_R_U": [point_projective_to_json(R_U) for R_U in self.vec_R_U],
            "x_final": field_to_json(self.x_final),
        }
        return json.dumps(dic)

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
    
    def to_bytes(self) -> bytes:
        return b''.join([
            g1_to_bytes(self.B_a),
            g1_to_bytes(self.B_t),
            g1_to_bytes(self.B_u),
            g1_list_to_bytes(self.vec_L_A),
            g1_list_to_bytes(self.vec_L_T),
            g1_list_to_bytes(self.vec_L_U),
            g1_list_to_bytes(self.vec_R_A),
            g1_list_to_bytes(self.vec_R_T),
            g1_list_to_bytes(self.vec_R_U),
            fr_to_bytes(self.x_final),
        ])
    
    @classmethod
    def from_bytes(cls: Type[T_SameMSMProof], b: BufReader, ell: int) -> T_SameMSMProof:
        log2_ell = int(log2(ell))
        return cls(
            B_a=b.read_g1(),
            B_t=b.read_g1(),
            B_u=b.read_g1(),
            vec_L_A=[b.read_g1() for _ in range(0, log2_ell)],
            vec_L_T=[b.read_g1() for _ in range(0, log2_ell)],
            vec_L_U=[b.read_g1() for _ in range(0, log2_ell)],
            vec_R_A=[b.read_g1() for _ in range(0, log2_ell)],
            vec_R_T=[b.read_g1() for _ in range(0, log2_ell)],
            vec_R_U=[b.read_g1() for _ in range(0, log2_ell)],
            x_final=b.read_fr(),
        )
