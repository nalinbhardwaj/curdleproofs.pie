import itertools
from crs import CurdleproofsCrs
from curdleproofs.same_perm import same_perm_verify
from curdleproofs.same_scalar import same_scalar_verify
from util import points_to_bytes, point_to_bytes
from msm_accumulator import MSMAccumulatorInefficient
from transcript import CurdleproofsTranscript
from typing import List, Tuple
from util import PointAffine, PointProjective, Fr, field_to_bytes
from py_ecc.optimized_bls12_381.optimized_curve import is_inf, add, neg, multiply

N_BLINDERS: int = 4

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


