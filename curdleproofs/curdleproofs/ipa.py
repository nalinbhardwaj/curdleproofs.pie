from math import log2
from curdleproofs.util import (
    field_from_json,
    field_to_json,
    point_projective_from_json,
    point_projective_to_json,
    points_projective_to_bytes,
    BufReader,
    g1_to_bytes,
    fr_to_bytes,
    g1_list_to_bytes,
    log2_int,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Tuple, Type, TypeVar
from curdleproofs.util import (
    field_to_bytes,
    invert,
    generate_blinders,
    inner_product,
    get_verification_scalars_bitstring,
)
from curdleproofs.msm_accumulator import MSMAccumulator, compute_MSM
from py_arkworks_bls12381 import G1Point, Scalar


def generate_ipa_blinders(c: List[Scalar], d: List[Scalar]) -> Tuple[List[Scalar], List[Scalar]]:
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

    assert inner_product(r, d) + inner_product(z, c) == Scalar(0)
    assert inner_product(r, z) == Scalar(0)

    return (r, z)


T_IPA = TypeVar("T_IPA", bound="IPA")


class IPA:
    def __init__(
        self,
        B_c: G1Point,
        B_d: G1Point,
        vec_L_C: List[G1Point],
        vec_R_C: List[G1Point],
        vec_L_D: List[G1Point],
        vec_R_D: List[G1Point],
        c_final: Scalar,
        d_final: Scalar,
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
        crs_G_vec: List[G1Point],
        crs_G_prime_vec: List[G1Point],
        crs_H: G1Point,
        C: G1Point,
        D: G1Point,
        z: Scalar,
        vec_c: List[Scalar],
        vec_d: List[Scalar],
        transcript: CurdleproofsTranscript,
    ) -> T_IPA:
        n = len(vec_c)
        lg_n = int(log2(n))
        if n != 2**lg_n:
            raise Exception("n != 2 ** lg_n, not a power of 2")
        if n != len(vec_d):
            raise Exception("len(vec_c) != len(vec_d)")

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
        H = crs_H * beta

        vec_L_C: List[G1Point] = []
        vec_R_C: List[G1Point] = []
        vec_L_D: List[G1Point] = []
        vec_R_D: List[G1Point] = []

        while len(vec_c) > 1:
            n //= 2
            # print("lens", len(vec_c), len(vec_d), len(crs_G_vec), len(crs_G_prime_vec))

            c_L, c_R = vec_c[:n], vec_c[n:]
            d_L, d_R = vec_d[:n], vec_d[n:]
            G_L, G_R = crs_G_vec[:n], crs_G_vec[n:]
            G_prime_L, G_prime_R = crs_G_prime_vec[:n], crs_G_prime_vec[n:]

            L_C = compute_MSM(G_R, c_L) + H * inner_product(c_L, d_R)
            L_D = compute_MSM(G_prime_L, d_R)
            R_C = compute_MSM(G_L, c_R) + H * inner_product(c_R, d_L)
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
                G_L[i] = G_L[i] + G_R[i] * gamma
                G_prime_L[i] = G_prime_L[i] + G_prime_R[i] * gamma_inv

            vec_c = c_L
            vec_d = d_L
            crs_G_vec = G_L
            crs_G_prime_vec = G_prime_L

        return cls(B_c, B_d, vec_L_C, vec_R_C, vec_L_D, vec_R_D, vec_c[0], vec_d[0])

    def verification_scalars(
        self, n: int, transcript: CurdleproofsTranscript
    ) -> Tuple[List[Scalar], List[Scalar], List[Scalar], List[Scalar]]:
        lg_n = len(self.vec_L_C)
        if lg_n >= 32:
            raise Exception("vec_L_C too large")
        elif n != 2**lg_n:
            raise Exception("n != 2 ** lg_n")

        verification_scalars_bitstring = get_verification_scalars_bitstring(n, lg_n)

        challenges: List[Scalar] = []
        for i in range(0, lg_n):
            transcript.append_list(
                b"ipa_loop",
                points_projective_to_bytes(
                    [self.vec_L_C[i], self.vec_L_D[i], self.vec_R_C[i], self.vec_R_D[i]]
                ),
            )
            challenges.append(transcript.get_and_append_challenge(b"ipa_gamma"))

        challenges_inv = [invert(c) for c in challenges]

        vec_s: List[Scalar] = []
        for i in range(0, n):
            vec_s.append(Scalar(1))
            for j in verification_scalars_bitstring[i]:
                vec_s[i] = vec_s[i] * challenges[j]

        vec_s_inv = [invert(s) for s in vec_s]

        return (challenges, challenges_inv, vec_s, vec_s_inv)

    def verify(
        self,
        crs_G_vec: List[G1Point],
        crs_H: G1Point,
        C: G1Point,
        D: G1Point,
        inner_prod: Scalar,
        vec_u: List[Scalar],
        transcript: CurdleproofsTranscript,
        msm_accumulator: MSMAccumulator,
    ):
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

        (vec_gamma, vec_gamma_inv, vec_s, vec_s_inv) = self.verification_scalars(
            n, transcript
        )

        vec_c_times_s = [self.c_final * s for s in vec_s]
        vec_rhs_scalars = vec_c_times_s + [self.c_final * self.d_final * beta]
        vec_G_H = crs_G_vec + [crs_H]

        H = crs_H * beta
        C_a = self.B_c + C * alpha + H * (alpha * alpha * inner_prod)

        point_lhs = compute_MSM(self.vec_L_C, vec_gamma) + C_a + compute_MSM(self.vec_R_C, vec_gamma_inv)

        msm_accumulator.accumulate_check(point_lhs, vec_G_H, vec_rhs_scalars)

        vec_d_div_s = [
            self.d_final * (s_inv_i * u_i) for (s_inv_i, u_i) in zip(vec_s_inv, vec_u)
        ]

        D_a = self.B_d + D * alpha
        point_lhs = compute_MSM(self.vec_L_D, vec_gamma) + D_a + compute_MSM(self.vec_R_D, vec_gamma_inv)
        msm_accumulator.accumulate_check(point_lhs, crs_G_vec, vec_d_div_s)

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
            c_final=field_from_json(json["c_final"]),
            d_final=field_from_json(json["d_final"]),
        )

    def to_bytes(self) -> bytes:
        return b''.join([
            g1_to_bytes(self.B_c),
            g1_to_bytes(self.B_d),
            g1_list_to_bytes(self.vec_L_C),
            g1_list_to_bytes(self.vec_R_C),
            g1_list_to_bytes(self.vec_L_D),
            g1_list_to_bytes(self.vec_R_D),
            fr_to_bytes(self.c_final),
            fr_to_bytes(self.d_final),
        ])

    @classmethod
    def from_bytes(cls: Type[T_IPA], b: BufReader, n: int) -> T_IPA:
        log2_n = log2_int(n)
        return cls(
            B_c=b.read_g1(),
            B_d=b.read_g1(),
            vec_L_C=[b.read_g1() for _ in range(0, log2_n)],
            vec_R_C=[b.read_g1() for _ in range(0, log2_n)],
            vec_L_D=[b.read_g1() for _ in range(0, log2_n)],
            vec_R_D=[b.read_g1() for _ in range(0, log2_n)],
            c_final=b.read_fr(),
            d_final=b.read_fr(),
        )
