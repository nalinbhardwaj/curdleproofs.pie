import itertools
from functools import reduce
from math import perm
import random

from crs import CurdleproofsCrs, get_random_point
from grand_prod import GrandProductProof
from ipa import generate_blinders
from transcript import CurdleproofsTranscript
from typing import List, Optional, Tuple, TypeVar
from util import PointAffine, PointProjective, Fr, field_to_bytes, points_projective_to_bytes
from msm_accumulator import MSMAccumulatorInefficient, SingleMSM
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, neg, Z1
from operator import mul as op_mul

T_SAME_PERM_PROOF = TypeVar('T_SAME_PERM_PROOF', bound="SamePermutationProof")

def get_permutation(vec_a: List[Fr], permutation: List[Fr]) -> List[Fr]:
  return [vec_a[i] for i in permutation]

class SamePermutationProof:
  def __init__(self,
    B: PointProjective,
    grand_prod_proof: GrandProductProof
  ) -> None:
    self.B = B
    self.grand_prod_proof = grand_prod_proof

  @classmethod
  def new(cls,
    crs_G_vec: List[PointProjective],
    crs_H_vec: List[PointProjective],
    crs_U: PointProjective,

    A: PointProjective,
    M: PointProjective,
    vec_a: List[Fr],

    permutation: List[int],
    vec_a_blinders: List[Fr],
    vec_m_blinders: List[Fr],

    transcript: CurdleproofsTranscript
  ) -> Tuple[Optional[T_SAME_PERM_PROOF], str]:
    n_blinders = len(vec_a_blinders)
    ell = len(crs_G_vec)

    transcript.append_list(b'same_perm_step1', points_projective_to_bytes([A, M]))
    transcript.append_list(b'same_perm_step1', list(map(field_to_bytes, vec_a)))
    alpha = transcript.get_and_append_challenge(b'same_perm_alpha')
    beta = transcript.get_and_append_challenge(b'same_perm_beta')

    vec_a_permuted = get_permutation(vec_a, permutation)
    permuted_polynomial_factors = [a + Fr(m) * alpha + beta for (a, m) in zip(vec_a_permuted, permutation)]
    gprod_result = reduce(op_mul, permuted_polynomial_factors, Fr.one())

    vec_beta_repeated = [beta] * ell
    B = add(add(A, multiply(M, int(alpha))), SingleMSM(crs_G_vec, list(map(int, vec_beta_repeated))).compute())

    vec_b_blinders = [vec_a_blinders[i] + alpha * vec_m_blinders[i] for i in range(0, n_blinders)]

    (grand_product_proof, err) = GrandProductProof.new(
      crs_G_vec=list(map(normalize, crs_G_vec)),
      crs_H_vec=list(map(normalize, crs_H_vec)),
      crs_U=crs_U,
      B=B,
      gprod_result=gprod_result,
      vec_b=permuted_polynomial_factors,
      vec_b_blinders=vec_b_blinders,
      transcript=transcript
    )

    if err is not None:
      return (None, err)

    return cls(B, grand_product_proof), ''

  def verify(self,
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
    msm_accumulator: MSMAccumulatorInefficient
  ) -> Tuple[bool, str]:
    ell = len(crs_G_vec)

    # Step 1
    transcript.append_list(b'same_perm_step1', points_projective_to_bytes([A, M]))
    transcript.append_list(b'same_perm_step1', list(map(field_to_bytes, vec_a)))
    
    alpha = transcript.get_and_append_challenge(b'same_perm_alpha')
    beta = transcript.get_and_append_challenge(b'same_perm_beta')

    # Step 2
    polynomial_factors = [a + Fr(i) * alpha + beta for (a, i) in zip(vec_a, range(0, ell))]
    gprod_result = reduce(op_mul, polynomial_factors, Fr.one())
    vec_beta_repeated = [beta] * ell
    msm_accumulator.accumulate_check(
      add(add(self.B, neg(A)), neg(multiply(M, int(alpha)))),
      crs_G_vec,
      list(map(int, vec_beta_repeated)),
    )

    (grand_prod_verify, err) = self.grand_prod_proof.verify(
      crs_G_vec=list(map(normalize, crs_G_vec)),
      crs_H_vec=list(map(normalize, crs_H_vec)),
      crs_U=crs_U,
      crs_G_sum=normalize(crs_G_sum),
      crs_H_sum=normalize(crs_H_sum),
      B=self.B,
      gprod_result=gprod_result,
      n_blinders=n_blinders,
      transcript=transcript,
      msm_accumulator=msm_accumulator
    )

    if not grand_prod_verify:
      return (False, err)
    
    return (True, '')

def test_same_permutation_proof():
  transcript_prover = CurdleproofsTranscript()

  n = 128
  n_blinders = 4
  ell = n - n_blinders

  crs_G_vec = [get_random_point() for _ in range(0, ell)]
  crs_H_vec = [get_random_point() for _ in range(0, n_blinders)]

  crs_U = get_random_point()
  crs_G_sum = reduce(add, crs_G_vec, Z1)
  crs_H_sum = reduce(add, crs_H_vec, Z1)

  vec_a_blinders = generate_blinders(n_blinders)
  vec_m_blinders = generate_blinders(n_blinders)

  permutation = list(range(0, ell))
  random.shuffle(permutation)

  vec_a = [Fr(random.randint(1, Fr.field_modulus - 1)) for _ in range(0, ell)]
  vec_a_permuted = get_permutation(vec_a, permutation)

  A = add(SingleMSM(crs_G_vec, list(map(int, vec_a_permuted))).compute(), SingleMSM(crs_H_vec, list(map(int, vec_a_blinders))).compute())
  M = add(SingleMSM(crs_G_vec, list(map(int, permutation))).compute(), SingleMSM(crs_H_vec, list(map(int, vec_m_blinders))).compute())

  (same_perm_proof, err) = SamePermutationProof.new(
    crs_G_vec=crs_G_vec,
    crs_H_vec=crs_H_vec,
    crs_U=crs_U,
    A=A,
    M=M,
    vec_a=vec_a,
    permutation=permutation,
    vec_a_blinders=vec_a_blinders,
    vec_m_blinders=vec_m_blinders,
    transcript=transcript_prover
  )

  print("Proof: ", same_perm_proof)
  print("Error: ", err)

  transcript_verifier = CurdleproofsTranscript()
  msm_accumulator = MSMAccumulatorInefficient()

  (verify, err) = same_perm_proof.verify(
    crs_G_vec=crs_G_vec,
    crs_H_vec=crs_H_vec,
    crs_U=crs_U,
    crs_G_sum=crs_G_sum,
    crs_H_sum=crs_H_sum,
    A=A,
    M=M,
    vec_a=vec_a,
    n_blinders=n_blinders,
    transcript=transcript_verifier,
    msm_accumulator=msm_accumulator
  )

  msm_verify = msm_accumulator.verify()

  print("Verify: ", verify)
  print("Error: ", err)
  print("MSM verify: ", msm_verify)

# test_same_permutation_proof()