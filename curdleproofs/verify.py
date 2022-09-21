import itertools
from crs import CurdleproofsCrs
from util import points_to_bytes, point_to_bytes
from msm_accumulator import MSMAccumulatorInefficient
from transcript import CurdleproofsTranscript
from typing import List, Tuple
from util import PointAffine, PointProjective, Fr, field_to_bytes
from py_ecc.optimized_bls12_381.optimized_curve import is_inf, add, neg, multiply
from operator import mul as op_mul

N_BLINDERS: int = 4

def ipa_verify(
  crs: CurdleproofsCrs,
  C: PointProjective,
  D: PointProjective,
  inner_prod: Fr,
  vec_u: List[Fr],
  B_c: PointProjective,
  B_d: PointProjective,
  vec_L_C: List[PointProjective],
  vec_R_C: List[PointProjective],
  vec_L_D: List[PointProjective],
  vec_R_D: List[PointProjective],
  c_final: Fr,
  d_final: Fr,
  transcript: CurdleproofsTranscript,
  msm_accumulator: MSMAccumulatorInefficient
):
  n = len(crs.vec_G)
  assert((n & (n-1) == 0) and n != 0, "n must be a power of 2")

  # Step 1
  transcript.append_list(b'ipa_step1', [C, D].map(points_to_bytes))
  transcript.append(b'ipa_step1', field_to_bytes(inner_prod))
  transcript.append_list(b'ipa_step1', [B_c, B_d].map(points_to_bytes))


  # transcript.append_list(b"ipa_step1", &[&C, &D]);
  #       transcript.append(b"ipa_step1", &z);
  #       transcript.append_list(b"ipa_step1", &[&self.B_c, &self.B_d]);
  #       let alpha = transcript.get_and_append_challenge(b"ipa_alpha");
  #       let beta = transcript.get_and_append_challenge(b"ipa_beta");


def grand_product_verify(crs: CurdleproofsCrs,
  B: PointProjective,
  C: PointProjective,
  r_p: Fr,
  gprod_result: Fr,
  n_blinders: int,
  transcript: CurdleproofsTranscript,
  msm_accumulator: MSMAccumulatorInefficient
) -> Tuple[bool, str]:
  ell = len(crs.vec_G)
  
  # Step 1
  transcript.append(b'gprod_step1', point_to_bytes(B))
  transcript.append(b'gprod_step1', field_to_bytes(gprod_result))
  alpha = transcript.get_and_append_challenge(b'gprod_alpha')

  # Step 2
  transcript.append(b'gprod_step1', C)
  transcript.append(b'gprod_step1', field_to_bytes(r_p))
  beta = transcript.get_and_append_challenge(b"gprod_beta");
  beta_inv = Fr.one() / beta

  # Step 3
  # Build `vec_u` for the optimization trick
  vec_u: List[Fr] = []
  pow_beta_inv = beta_inv
  for _ in range(0, ell):
    vec_u.append(pow_beta_inv)
    pow_beta_inv *= beta_inv
    
  vec_u.extend([beta_inv ** (ell + 1) for _ in range(0, n_blinders)])

  # Compute D
  D = add(add(B, neg(multiply(crs.G_sum, beta_inv))), multiply(crs.H_sum, alpha))

  # Step 4
  # Build G
  vec_G = crs.vec_G.copy()
  vec_G.extend(crs.vec_H)

  inner_prod = r_p * (beta ** (ell + 1)) + gprod_result * (beta ** ell) - Fr.one()

  return ipa_verify(
    crs,
    C,
    D,
    inner_prod,
    vec_u,
    transcript,
    msm_accumulator
  )

def same_perm_verify(crs: CurdleproofsCrs,
  A: PointProjective,
  B: PointProjective,
  M: PointProjective,
  vec_a: List[Fr],
  n_blinders: int,
  transcript: CurdleproofsTranscript,
  msm_accumulator: MSMAccumulatorInefficient
) -> Tuple[bool, str]:
  ell = len(crs.vec_G)

  # Step 1
  transcript.append_list(b'same_perm_step1', [A, M].map(points_to_bytes))
  transcript.append_list(b'same_perm_step1', [vec_a].map(points_to_bytes))
  
  alpha = transcript.get_and_append_challenge(b'same_perm_alpha')
  beta = transcript.get_and_append_challenge(b'same_perm_beta')

  # Step 2
  range_as_fr = [Fr(i) for i in range(0, ell)]
  polynomial_factors = [a + i * alpha + beta for (a, i) in zip(vec_a, range_as_fr)]
  gprod_result = itertools.accumulate(polynomial_factors, op_mul)
  vec_beta_repeated = [beta] * ell
  msm_accumulator.accumulate_check(
    add(add(B, neg(A)), neg(multiply(M, alpha))),
    vec_beta_repeated,
    crs.vec_G
  )

  return grand_product_verify(
    crs,
    B,
    gprod_result,
    n_blinders,
    transcript,
    msm_accumulator
  )

def same_scalar_verify():
  pass


def verify(crs: CurdleproofsCrs,
  vec_R: List[PointAffine],
  vec_S: List[PointAffine],
  vec_T: List[PointAffine],
  vec_U: List[PointAffine],
  M: PointProjective,
  A: PointProjective,
  B: PointProjective,
) -> Tuple[bool, str]:
  ell = len(vec_R)

  transcript = CurdleproofsTranscript()

  msm_accumulator = MSMAccumulatorInefficient()

  if is_inf(vec_T[0]):
    return (False, 'randomizer was the zero element')

  # Step 1
  transcript.append_list(b'curdleproofs_step1', [vec_R, vec_S, vec_T, vec_U].map(points_to_bytes))
  transcript.append(b'curdleproofs_step1', point_to_bytes(M))
  vec_a = transcript.get_and_append_challenges(b'curdleproofs_vec_a', ell)

  # Step 2
  # Verify the grand product proof
  same_perm_verify(
    crs,
    A,
    B,
    M,
    vec_a,
    N_BLINDERS,
    transcript,
    msm_accumulator
  );
  
  # Step 3
  same_scalar_verify()

  # Step 4


