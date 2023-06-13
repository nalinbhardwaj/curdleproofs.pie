import json
from math import log2
import random
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
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Optional, Tuple, Type, TypeVar
from curdleproofs.util import (
    PointProjective,
    Fr,
    field_to_bytes,
    invert,
    generate_blinders,
    inner_product,
    get_verification_scalars_bitstring,
)
from curdleproofs.msm_accumulator import MSMAccumulator, compute_MSM
from py_ecc.optimized_bls12_381.optimized_curve import (
    curve_order,
    G1,
    multiply,
    normalize,
    add,
    neg,
)


def generate_ipa_blinders(c: List[Fr], d: List[Fr]) -> Tuple[List[Fr], List[Fr]]:
    n = len(c)

    r = generate_blinders(n)
    z = generate_blinders(n - 2)

    omega = inner_product(r, d) + inner_product(z[: n - 2], c[: n - 2])
    delta = inner_product(r[: n - 2], z[: n - 2])

    inv_c = invert(c[n - 2])

    last_z = (r[n - 2] * inv_c * omega - delta) * invert(
        -r[n - 2] * inv_c * c[n - 1] + r[n - 1]
    )
    penultimate_z = -inv_c * (last_z * c[n - 1] + omega)

    z += [penultimate_z, last_z]

    assert inner_product(r, d) + inner_product(z, c) == Fr.zero()
    assert inner_product(r, z) == Fr.zero()

    return (r, z)


T_IPA = TypeVar("T_IPA", bound="IPA")


class IPA:
    def __init__(
        self,
        B_c: PointProjective,
        B_d: PointProjective,
        vec_L_C: List[PointProjective],
        vec_R_C: List[PointProjective],
        vec_L_D: List[PointProjective],
        vec_R_D: List[PointProjective],
        c_final: Fr,
        d_final: Fr,
    ) -> None:
        self.B_c = B_c
        self.B_d = B_d
        self.vec_L_C = vec_L_C
        self.vec_R_C = vec_R_C
        self.vec_L_D = vec_L_D
        self.vec_R_D = vec_R_D
        self.c_final = c_final
        self.d_final = d_final

    @classmethod
    def new(
        cls: Type[T_IPA],
        crs_G_vec: List[PointProjective],
        crs_G_prime_vec: List[PointProjective],
        crs_H: PointProjective,
        C: PointProjective,
        D: PointProjective,
        z: Fr,
        vec_c: List[Fr],
        vec_d: List[Fr],
        transcript: CurdleproofsTranscript,
    ) -> Tuple[Optional[T_IPA], Optional[str]]:
        n = len(vec_c)
        lg_n = int(log2(n))
        if n != 2**lg_n:
            return (None, "n != 2 ** lg_n, not a power of 2")
        if n != len(vec_d):
            return (None, "len(vec_c) != len(vec_d)")

        (vec_r_c, vec_r_d) = generate_ipa_blinders(vec_c, vec_d)

        B_c = compute_MSM(crs_G_vec, vec_r_c)
        B_d = compute_MSM(crs_G_prime_vec, vec_r_d)

        transcript.append_list(b"ipa_step1", points_projective_to_bytes([C, D]))
        transcript.append(b"ipa_step1", field_to_bytes(z))
        transcript.append_list(b"ipa_step1", points_projective_to_bytes([B_c, B_d]))

        alpha = transcript.get_and_append_challenge(b"ipa_alpha")
        beta = transcript.get_and_append_challenge(b"ipa_beta")

        for i in range(0, n):
            vec_c[i] = vec_r_c[i] + alpha * vec_c[i]
            vec_d[i] = vec_r_d[i] + alpha * vec_d[i]
        H = multiply(crs_H, int(beta))

        vec_L_C: List[PointProjective] = []
        vec_R_C: List[PointProjective] = []
        vec_L_D: List[PointProjective] = []
        vec_R_D: List[PointProjective] = []

        while len(vec_c) > 1:
            n //= 2
            # print("lens", len(vec_c), len(vec_d), len(crs_G_vec), len(crs_G_prime_vec))

            c_L, c_R = vec_c[:n], vec_c[n:]
            d_L, d_R = vec_d[:n], vec_d[n:]
            G_L, G_R = crs_G_vec[:n], crs_G_vec[n:]
            G_prime_L, G_prime_R = crs_G_prime_vec[:n], crs_G_prime_vec[n:]

            L_C = add(compute_MSM(G_R, c_L), multiply(H, int(inner_product(c_L, d_R))))
            L_D = compute_MSM(G_prime_L, d_R)
            R_C = add(compute_MSM(G_L, c_R), multiply(H, int(inner_product(c_R, d_L))))
            R_D = compute_MSM(G_prime_R, d_L)

            vec_L_C.append(L_C)
            vec_R_C.append(R_C)
            vec_L_D.append(L_D)
            vec_R_D.append(R_D)

            transcript.append_list(
                b"ipa_loop", points_projective_to_bytes([L_C, L_D, R_C, R_D])
            )
            gamma = transcript.get_and_append_challenge(b"ipa_gamma")
            gamma_inv = invert(gamma)

            for i in range(0, n):
                c_L[i] += gamma_inv * c_R[i]
                d_L[i] += gamma * d_R[i]
                G_L[i] = add(G_L[i], multiply(G_R[i], int(gamma)))
                G_prime_L[i] = add(G_prime_L[i], multiply(G_prime_R[i], int(gamma_inv)))

            vec_c = c_L
            vec_d = d_L
            crs_G_vec = G_L
            crs_G_prime_vec = G_prime_L

        return (
            cls(B_c, B_d, vec_L_C, vec_R_C, vec_L_D, vec_R_D, vec_c[0], vec_d[0]),
            None,
        )

    def verification_scalars(
        self, n: int, transcript: CurdleproofsTranscript
    ) -> Tuple[Tuple[List[Fr], List[Fr], List[Fr], List[Fr]], Optional[str]]:
        lg_n = len(self.vec_L_C)
        if lg_n >= 32:
            return (([], [], [], []), "vec_L_C too large")
        elif n != 2**lg_n:
            return (([], [], [], []), "n != 2 ** lg_n")

        verification_scalars_bitstring = get_verification_scalars_bitstring(n, lg_n)

        challenges: List[Fr] = []
        for i in range(0, lg_n):
            transcript.append_list(
                b"ipa_loop",
                points_projective_to_bytes(
                    [self.vec_L_C[i], self.vec_L_D[i], self.vec_R_C[i], self.vec_R_D[i]]
                ),
            )
            challenges.append(transcript.get_and_append_challenge(b"ipa_gamma"))

        challenges_inv = [invert(c) for c in challenges]

        vec_s: List[Fr] = []
        for i in range(0, n):
            vec_s.append(Fr.one())
            for j in verification_scalars_bitstring[i]:
                vec_s[i] = vec_s[i] * challenges[j]

        vec_s_inv = [invert(s) for s in vec_s]

        return ((challenges, challenges_inv, vec_s, vec_s_inv), None)

    def verify(
        self,
        crs_G_vec: List[PointProjective],
        crs_H: PointProjective,
        C: PointProjective,
        D: PointProjective,
        inner_prod: Fr,
        vec_u: List[Fr],
        transcript: CurdleproofsTranscript,
        msm_accumulator: MSMAccumulator,
    ) -> Tuple[bool, str]:
        n = len(crs_G_vec)
        # assert(((n != 0) and (n & (n-1) == 0)), "n must be a power of 2")

        # Step 1
        transcript.append_list(b"ipa_step1", points_projective_to_bytes([C, D]))
        transcript.append(b"ipa_step1", field_to_bytes(inner_prod))
        transcript.append_list(
            b"ipa_step1", points_projective_to_bytes([self.B_c, self.B_d])
        )

        alpha = transcript.get_and_append_challenge(b"ipa_alpha")
        beta = transcript.get_and_append_challenge(b"ipa_beta")

        ((vec_gamma, vec_gamma_inv, vec_s, vec_s_inv), err) = self.verification_scalars(
            n, transcript
        )
        if err is not None:
            return (False, err)

        vec_c_times_s = [self.c_final * s for s in vec_s]
        vec_rhs_scalars = vec_c_times_s + [self.c_final * self.d_final * beta]
        vec_G_H = crs_G_vec + [crs_H]

        H = multiply(crs_H, int(beta))
        C_a: PointProjective = add(
            add(self.B_c, multiply(C, int(alpha))),
            multiply(H, int(alpha * alpha * inner_prod)),
        )

        point_lhs = add(
            add(compute_MSM(self.vec_L_C, vec_gamma), C_a),
            compute_MSM(self.vec_R_C, vec_gamma_inv),
        )

        msm_accumulator.accumulate_check(point_lhs, vec_G_H, vec_rhs_scalars)

        vec_d_div_s = [
            self.d_final * (s_inv_i * u_i) for (s_inv_i, u_i) in zip(vec_s_inv, vec_u)
        ]

        D_a = add(self.B_d, multiply(D, int(alpha)))
        point_lhs = add(
            add(compute_MSM(self.vec_L_D, vec_gamma), D_a),
            compute_MSM(self.vec_R_D, vec_gamma_inv),
        )
        msm_accumulator.accumulate_check(point_lhs, crs_G_vec, vec_d_div_s)

        return (True, "")

    def to_json(self):
        return {
            "B_c": point_projective_to_json(self.B_c),
            "B_d": point_projective_to_json(self.B_d),
            "vec_L_C": [point_projective_to_json(p) for p in self.vec_L_C],
            "vec_R_C": [point_projective_to_json(p) for p in self.vec_R_C],
            "vec_L_D": [point_projective_to_json(p) for p in self.vec_L_D],
            "vec_R_D": [point_projective_to_json(p) for p in self.vec_R_D],
            "c_final": field_to_json(self.c_final),
            "d_final": field_to_json(self.d_final),
        }

    @classmethod
    def from_json(cls: Type[T_IPA], json) -> T_IPA:
        return cls(
            B_c=point_projective_from_json(json["B_c"]),
            B_d=point_projective_from_json(json["B_d"]),
            vec_L_C=[point_projective_from_json(p) for p in json["vec_L_C"]],
            vec_R_C=[point_projective_from_json(p) for p in json["vec_R_C"]],
            vec_L_D=[point_projective_from_json(p) for p in json["vec_L_D"]],
            vec_R_D=[point_projective_from_json(p) for p in json["vec_R_D"]],
            c_final=field_from_json(json["c_final"], Fr),
            d_final=field_from_json(json["d_final"], Fr),
        )
