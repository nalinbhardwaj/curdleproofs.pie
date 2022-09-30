from crs import CurdleproofsCrs
from ipa import ipa_verify
from util import points_to_bytes, point_to_bytes
from transcript import CurdleproofsTranscript
from typing import List, Tuple
from util import PointAffine, PointProjective, Fr, field_to_bytes
from msm_accumulator import MSMAccumulatorInefficient
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, neg

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
