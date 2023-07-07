from functools import reduce
import json
import operator
import random
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.ipa import IPA, generate_blinders, inner_product
from curdleproofs.util import (
    field_from_json,
    field_to_json,
    invert,
    point_projective_from_json,
    point_projective_to_bytes,
    get_random_point,
    point_projective_to_json,
    BufReader,
    g1_to_bytes,
    fr_to_bytes,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Optional, Tuple, TypeVar, Type
from curdleproofs.util import (
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

    @classmethod
    def new(
        cls: Type[T_GrandProductProof],
        crs_G_vec: List[PointProjective],
        crs_H_vec: List[PointProjective],
        crs_U: PointProjective,
        B: PointProjective,
        gprod_result: Fr,
        vec_b: List[Fr],
        vec_b_blinders: List[Fr],
        transcript: CurdleproofsTranscript,
    ) -> Tuple[Optional[T_GrandProductProof], Optional[str]]:
        n_blinders = len(vec_b_blinders)
        ell = len(crs_G_vec)
        n = ell + n_blinders

        transcript.append(b"gprod_step1", point_projective_to_bytes(B))
        transcript.append(b"gprod_step1", field_to_bytes(gprod_result))
        alpha = transcript.get_and_append_challenge(b"gprod_alpha")

        # Step 2
        vec_c = [Fr.one()]
        for i in range(0, ell - 1):
            vec_c.append(vec_c[i] * vec_b[i])

        vec_c_blinders = generate_blinders(n_blinders)
        C = add(compute_MSM(crs_G_vec, vec_c), compute_MSM(crs_H_vec, vec_c_blinders))

        vec_r_b_plus_alpha = [r_b_i + alpha for r_b_i in vec_b_blinders]
        r_p = inner_product(vec_r_b_plus_alpha, vec_c_blinders)

        transcript.append(b"gprod_step2", point_projective_to_bytes(C))
        transcript.append(b"gprod_step2", field_to_bytes(r_p))
        beta = transcript.get_and_append_challenge(b"gprod_beta")
        beta_inv = invert(beta)

        pow_beta_inv = beta_inv
        vec_G_prime: List[PointProjective] = []
        for G_i in crs_G_vec:
            G_prime = multiply(G_i, int(pow_beta_inv))
            vec_G_prime.append(G_prime)
            pow_beta_inv *= beta_inv

        vec_H_prime = [multiply(H_i, int(beta_inv ** (ell + 1))) for H_i in crs_H_vec]

        vec_b_prime: List[Fr] = []
        pow_beta = beta
        for b_i in vec_b:
            vec_b_prime.append(b_i * pow_beta)
            pow_beta *= beta

        vec_d: List[Fr] = []
        pow_beta = Fr.one()
        vec_beta_powers: List[Fr] = []
        for b_prime_i in vec_b_prime:
            vec_d.append(b_prime_i - pow_beta)
            vec_beta_powers.append(pow_beta)
            pow_beta *= beta

        vec_d_blinders = [(beta ** (ell + 1)) * r_b_i for r_b_i in vec_r_b_plus_alpha]

        vec_alphabeta = [alpha * (beta ** (ell + 1)) for _ in range(n_blinders)]
        D = add(
            add(B, neg(compute_MSM(vec_G_prime, vec_beta_powers))),
            compute_MSM(vec_H_prime, vec_alphabeta),
        )

        vec_G = crs_G_vec + crs_H_vec
        vec_G_prime += vec_H_prime

        inner_prod = r_p * (beta ** (ell + 1)) + gprod_result * (beta**ell) - Fr.one()

        vec_c += vec_c_blinders
        vec_d += vec_d_blinders

        # print("inner_prod", inner_prod)
        # print("computed", inner_product(vec_c, vec_d))

        assert inner_product(vec_c, vec_d) == inner_prod
        assert eq(compute_MSM(vec_G, vec_c), C)
        assert eq(compute_MSM(vec_G_prime, vec_d), D)

        (ipa_proof, err) = IPA.new(
            crs_G_vec=vec_G,
            crs_G_prime_vec=vec_G_prime,
            crs_H=crs_U,
            C=C,
            D=D,
            z=inner_prod,
            vec_c=vec_c,
            vec_d=vec_d,
            transcript=transcript,
        )

        if ipa_proof is None:
            return None, err

        return cls(C, r_p, ipa_proof), None

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

    def to_json(self):
        return {
            "C": point_projective_to_json(self.C),
            "r_p": field_to_json(self.r_p),
            "ipa_proof": self.ipa_proof.to_json(),
        }

    @classmethod
    def from_json(cls: Type[T_GrandProductProof], json) -> T_GrandProductProof:
        return cls(
            C=point_projective_from_json(json["C"]),
            r_p=field_from_json(json["r_p"], Fr),
            ipa_proof=IPA.from_json(json["ipa_proof"]),
        )
    
    def to_bytes(self) -> bytes:
        return b''.join([
            g1_to_bytes(self.C),
            fr_to_bytes(self.r_p),
            self.ipa_proof.to_bytes(),
        ])

    @classmethod
    def from_bytes(cls: Type[T_GrandProductProof], b: BufReader, n: int) -> T_GrandProductProof:
        return cls(
            C=b.read_g1(),
            r_p=b.read_fr(),
            ipa_proof=IPA.from_bytes(b, n),
        )
