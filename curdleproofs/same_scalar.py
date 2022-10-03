import itertools
import random

from crs import CurdleproofsCrs, get_random_point
from commitment import GroupCommitment
from util import points_projective_to_bytes
from transcript import CurdleproofsTranscript
from typing import List, Tuple, TypeVar
from util import PointAffine, PointProjective, Fr, field_to_bytes
from msm_accumulator import MSMAccumulatorInefficient
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, neg
from operator import mul as op_mul

T_SameScalarProof = TypeVar('T_SameScalarProof', bound="SameScalarProof")
class SameScalarProof:
  def __init__(self, cm_A: GroupCommitment, cm_B: GroupCommitment, z_k: Fr, z_t: Fr, z_u: Fr) -> None:
    self.cm_A = cm_A
    self.cm_B = cm_B
    self.z_k = z_k
    self.z_t = z_t
    self.z_u = z_u

  @classmethod
  def new(cls,
    crs_G_t: PointProjective,
    crs_G_u: PointProjective,
    crs_H: PointProjective,

    R: PointProjective,
    S: PointProjective,
    cm_T: GroupCommitment,
    cm_U: GroupCommitment,
    k: Fr,
    r_t: Fr,
    r_u: Fr,
    transcript: CurdleproofsTranscript
  ) -> T_SameScalarProof:
    r_a = Fr(random.randint(1, Fr.field_modulus))
    r_b = Fr(random.randint(1, Fr.field_modulus))
    r_k = Fr(random.randint(1, Fr.field_modulus))

    cm_A = GroupCommitment.new(crs_G_t, crs_H, multiply(R, int(r_k)), r_a)
    cm_B = GroupCommitment.new(crs_G_u, crs_H, multiply(S, int(r_k)), r_b)

    transcript.append_list(b'sameexp_points', points_projective_to_bytes([R, S, cm_T.T_1, cm_T.T_2, cm_U.T_1, cm_U.T_2, cm_A.T_1, cm_A.T_2, cm_B.T_1, cm_B.T_2]))
    alpha = transcript.get_and_append_challenge(b'same_scalar_alpha')

    z_k = r_k + k * alpha
    z_t = r_a + r_t * alpha
    z_u = r_b + r_u * alpha

    return cls(cm_A, cm_B, z_k, z_t, z_u)

  def verify(self,
    crs_G_t: PointProjective,
    crs_G_u: PointProjective,
    crs_H: PointProjective,

    R: PointProjective,
    S: PointProjective,
    cm_T: GroupCommitment,
    cm_U: GroupCommitment,
    transcript: CurdleproofsTranscript
  ) -> Tuple[bool, str]:
    transcript.append_list(b'sameexp_points', points_projective_to_bytes([R, S, cm_T.T_1, cm_T.T_2, cm_U.T_1, cm_U.T_2, self.cm_A.T_1, self.cm_A.T_2, self.cm_B.T_1, self.cm_B.T_2]))

    alpha = transcript.get_and_append_challenge(b'same_scalar_alpha');
    expected_1 = GroupCommitment.new(crs_G_t, crs_H, multiply(R, int(self.z_k)), self.z_t);
    expected_2 = GroupCommitment.new(crs_G_u, crs_H, multiply(S, int(self.z_k)), self.z_u);

    computed_1 = self.cm_A + (cm_T * alpha)
    computed_2 = self.cm_B + (cm_U * alpha)

    if expected_1 == computed_1 and expected_2 == computed_2:
      return (True, "")
    else:
      return (False, "Failure")

def test_same_scalar_arg():
  transcript_prover = CurdleproofsTranscript()

  crs_G_t = get_random_point()
  crs_G_u = get_random_point()
  crs_H = get_random_point()

  R = get_random_point()
  S = get_random_point()

  k = Fr(random.randint(1, Fr.field_modulus))
  r_t = Fr(random.randint(1, Fr.field_modulus))
  r_u = Fr(random.randint(1, Fr.field_modulus))

  cm_T = GroupCommitment.new(crs_G_t, crs_H, multiply(R, int(k)), r_t)
  cm_U = GroupCommitment.new(crs_G_u, crs_H, multiply(S, int(k)), r_u)

  proof = SameScalarProof.new(
    crs_G_t=crs_G_t,
    crs_G_u=crs_G_u,
    crs_H=crs_H,
    R=R,
    S=S,
    cm_T=cm_T,
    cm_U=cm_U,
    k=k,
    r_t=r_t,
    r_u=r_u,
    transcript=transcript_prover
  )

  print("proof", proof)

  transcript_verifier = CurdleproofsTranscript()
  (res, err) = proof.verify(
    crs_G_t=crs_G_t,
    crs_G_u=crs_G_u,
    crs_H=crs_H,
    R=R,
    S=S,
    cm_T=cm_T,
    cm_U=cm_U,
    transcript=transcript_verifier
  )
  print("res", res)
  print("err", err)

# test_same_scalar_arg()