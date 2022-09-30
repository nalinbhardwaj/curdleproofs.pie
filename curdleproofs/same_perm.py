import itertools

from crs import CurdleproofsCrs
from grand_prod import grand_product_verify
from ipa import ipa_verify
from util import points_to_bytes, point_to_bytes
from transcript import CurdleproofsTranscript
from typing import List, Tuple
from util import PointAffine, PointProjective, Fr, field_to_bytes
from msm_accumulator import MSMAccumulatorInefficient
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, neg
from operator import mul as op_mul

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