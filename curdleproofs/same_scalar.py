import itertools

from crs import CurdleproofsCrs
from curdleproofs.commitment import GroupCommitment
from grand_prod import grand_product_verify
from ipa import ipa_verify
from util import points_to_bytes, point_to_bytes
from transcript import CurdleproofsTranscript
from typing import List, Tuple
from util import PointAffine, PointProjective, Fr, field_to_bytes
from msm_accumulator import MSMAccumulatorInefficient
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, neg
from operator import mul as op_mul

class SameScalarProof:
  def __init__(self, cm_A: GroupCommitment, cm_B: GroupCommitment, z_k: Fr, z_t: Fr, z_u: Fr) -> None:
    self.cm_A = cm_A
    self.cm_B = cm_B
    self.z_k = z_k
    self.z_t = z_t
    self.z_u = z_u

  def verify(self,
    crs: CurdleproofsCrs,
    R: PointProjective,
    S: PointProjective,
    cm_T: GroupCommitment,
    cm_U: GroupCommitment,
    transcript: CurdleproofsTranscript) -> Tuple[bool, str]:
    transcript.append_list(b'sameexp_points', [R, S, cm_T.T_1, cm_T.T_2, cm_U.T_1, cm_U.T_2, self.cm_A.T_1, self.cm_A.T_2, self.cm_B.T_1, self.cm_B.T_2].map(points_to_bytes))

    alpha = transcript.get_and_append_challenge(b'same_scalar_alpha');
    expected_1 = GroupCommitment(crs.G_t, crs.H, multiply(R, int(self.z_k)), self.z_t);
    expected_2 = GroupCommitment(crs.G_t, crs.H, multiply(S, int(self.z_k)), self.z_u);

    computed_1 = self.cm_A + cm_T * alpha
    computed_2 = self.cm_B + cm_U * alpha

    if expected_1.T_1 == computed_1.T_1 and expected_1.T_2 == computed_1.T_2 and expected_2.T_1 == computed_2.T_1 and expected_2.T_2 == computed_2.T_2:
      return (True, "")
    else:
      return (False, "Failure")