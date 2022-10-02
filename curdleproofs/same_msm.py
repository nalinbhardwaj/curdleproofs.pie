from math import log2
import random
from crs import CurdleproofsCrs, get_random_point
from ipa import generate_blinders, get_verification_scalars_bitstring
from msm_accumulator import SingleMSM
from util import affine_to_projective, point_affine_to_bytes, point_projective_to_bytes, points_affine_to_bytes, points_projective_to_bytes
from transcript import CurdleproofsTranscript
from typing import List, Optional, Tuple, Type, TypeVar
from util import PointAffine, PointProjective, Fr, field_to_bytes, invert
from msm_accumulator import MSMAccumulatorInefficient
from py_ecc.optimized_bls12_381.optimized_curve import curve_order, G1, multiply, normalize, add, neg

T_SameMSMProof = TypeVar('T_SameMSMProof', bound="SameMSMProof")

class SameMSMProof:
  def __init__(self,
    B_a: PointProjective,
    B_t: PointProjective,
    B_u: PointProjective,

    vec_L_A: List[PointProjective],
    vec_L_T: List[PointProjective],
    vec_L_U: List[PointProjective],
    vec_R_A: List[PointProjective],
    vec_R_T: List[PointProjective],
    vec_R_U: List[PointProjective],

    x_final: Fr
  ) -> None:
    self.B_a = B_a
    self.B_t = B_t
    self.B_u = B_u
    self.vec_L_A = vec_L_A
    self.vec_L_T = vec_L_T
    self.vec_L_U = vec_L_U
    self.vec_R_A = vec_R_A
    self.vec_R_T = vec_R_T
    self.vec_R_U = vec_R_U
    self.x_final = x_final
  
  @classmethod
  def new(cls,
    crs_G_vec: List[PointAffine],

    A: PointProjective,
    Z_t: PointProjective,
    Z_u: PointProjective,
    vec_T: List[PointAffine],
    vec_U: List[PointAffine],

    vec_x: List[Fr],

    transcript: CurdleproofsTranscript,
  ) -> T_SameMSMProof:
    n = len(vec_x)
    lg_n = int(log2(n))
    assert 2**lg_n == n

    vec_L_T: List[PointProjective] = []
    vec_R_T: List[PointProjective] = []
    vec_L_U: List[PointProjective] = []
    vec_R_U: List[PointProjective] = []
    vec_L_A: List[PointProjective] = []
    vec_R_A: List[PointProjective] = []

    vec_r = generate_blinders(n)

    B_a = SingleMSM(list(map(affine_to_projective, crs_G_vec)), list(map(int, vec_r))).compute()
    B_t = SingleMSM(list(map(affine_to_projective, vec_T)), list(map(int, vec_r))).compute()
    B_u = SingleMSM(list(map(affine_to_projective, vec_U)), list(map(int, vec_r))).compute()

    transcript.append_list(b'same_msm_step1', points_projective_to_bytes([A, Z_t, Z_u]))
    transcript.append_list(b'same_msm_step1', points_affine_to_bytes(vec_T + vec_U))
    transcript.append_list(b'same_msm_step1', points_projective_to_bytes([B_a, B_t, B_u]))
    alpha = transcript.get_and_append_challenge(b'same_msm_alpha')

    for i in range(0, n):
      vec_x[i] = vec_r[i] + (alpha * vec_x[i])
    
    while len(vec_x) > 1:
      n //= 2

      x_L, x_R = vec_x[:n], vec_x[n:]
      T_L, T_R = vec_T[:n], vec_T[n:]
      U_L, U_R = vec_U[:n], vec_U[n:]
      G_L, G_R = crs_G_vec[:n], crs_G_vec[n:]

      L_A = SingleMSM(list(map(affine_to_projective, G_R)), list(map(int, x_L))).compute()
      L_T = SingleMSM(list(map(affine_to_projective, T_R)), list(map(int, x_L))).compute()
      L_U = SingleMSM(list(map(affine_to_projective, U_R)), list(map(int, x_L))).compute()
      R_A = SingleMSM(list(map(affine_to_projective, G_L)), list(map(int, x_R))).compute()
      R_T = SingleMSM(list(map(affine_to_projective, T_L)), list(map(int, x_R))).compute()
      R_U = SingleMSM(list(map(affine_to_projective, U_L)), list(map(int, x_R))).compute()

      vec_L_A.append(L_A)
      vec_L_T.append(L_T)
      vec_L_U.append(L_U)
      vec_R_A.append(R_A)
      vec_R_T.append(R_T)
      vec_R_U.append(R_U)

      transcript.append_list(b'same_msm_loop', points_projective_to_bytes([L_A, L_T, L_U, R_A, R_T, R_U]))
      gamma = transcript.get_and_append_challenge(b'same_msm_gamma')
      gamma_inv = invert(gamma)

      for i in range(0, n):
        x_L[i] += gamma_inv * x_R[i]
        T_L[i] = normalize(add(affine_to_projective(T_L[i]), multiply(affine_to_projective(T_R[i]), int(gamma))))
        U_L[i] = normalize(add(affine_to_projective(U_L[i]), multiply(affine_to_projective(U_R[i]), int(gamma))))
        G_L[i] = normalize(add(affine_to_projective(G_L[i]), multiply(affine_to_projective(G_R[i]), int(gamma))))
      
      vec_x = x_L
      vec_T = T_L
      vec_U = U_L
      crs_G_vec = G_L

    return cls(
      B_a=B_a,
      B_t=B_t,
      B_u=B_u,
      vec_L_A=vec_L_A,
      vec_L_T=vec_L_T,
      vec_L_U=vec_L_U,
      vec_R_A=vec_R_A,
      vec_R_T=vec_R_T,
      vec_R_U=vec_R_U,
      x_final=vec_x[0]
    )

  def verification_scalars(self, n: int, transcript: CurdleproofsTranscript) -> Tuple[Optional[Tuple[List[Fr], List[Fr], List[Fr]]], str]:
    lg_n = len(self.vec_L_A)
    if lg_n >= 32:
      return None, 'lg_n >= 32'
    if 2**lg_n != n:
      return None, '2**lg_n != n'
    
    bitstring = get_verification_scalars_bitstring(n, lg_n)

    challenges: List[Fr] = []
    for i in range(0, lg_n):
      transcript.append_list(b'same_msm_loop', points_projective_to_bytes([self.vec_L_A[i], self.vec_L_T[i], self.vec_L_U[i], self.vec_R_A[i], self.vec_R_T[i], self.vec_R_U[i]]))
      challenges.append(transcript.get_and_append_challenge(b'same_msm_gamma'))

    challenges_inv = list(map(invert, challenges))

    vec_s: List[Fr] = []
    for i in range(0, n):
      vec_s.append(Fr.one())
      for j in bitstring[i]:
        vec_s[i] *= challenges[j]

    return (challenges, challenges_inv, vec_s), ''
  
  def verify(self,
    crs_G_vec: List[PointProjective],

    A: PointProjective,
    Z_t: PointProjective,
    Z_u: PointProjective,
    vec_T: List[PointAffine],
    vec_U: List[PointAffine],
    transcript: CurdleproofsTranscript,
    msm_accumulator: MSMAccumulatorInefficient
  ) -> Tuple[bool, str]:
    n = len(vec_T)

    transcript.append_list(b'same_msm_step1', points_projective_to_bytes([A, Z_t, Z_u]))
    transcript.append_list(b'same_msm_step1', points_affine_to_bytes(vec_T + vec_U))
    transcript.append_list(b'same_msm_step1', points_projective_to_bytes([self.B_a, self.B_t, self.B_u]))
    alpha = transcript.get_and_append_challenge(b'same_msm_alpha')

    ((vec_gamma, vec_gamma_inv, vec_s), err) = self.verification_scalars(n, transcript)

    if vec_gamma is None:
      return False, err
    
    vec_x_times_s = [self.x_final * s_i for s_i in vec_s]

    A_a = add(self.B_a, multiply(A, int(alpha)))
    Z_t_a = add(self.B_t, multiply(Z_t, int(alpha)))
    Z_u_a = add(self.B_u, multiply(Z_u, int(alpha)))

    point_lhs = add(add(SingleMSM(self.vec_L_A, list(map(int, vec_gamma))).compute(), A_a), SingleMSM(self.vec_R_A, list(map(int, vec_gamma_inv))).compute())
    msm_accumulator.accumulate_check(point_lhs, crs_G_vec, list(map(int, vec_x_times_s)))

    point_lhs = add(add(SingleMSM(self.vec_L_T, list(map(int, vec_gamma))).compute(), Z_t_a), SingleMSM(self.vec_R_T, list(map(int, vec_gamma_inv))).compute())
    msm_accumulator.accumulate_check(point_lhs, list(map(affine_to_projective, vec_T)), list(map(int, vec_x_times_s)))

    point_lhs = add(add(SingleMSM(self.vec_L_U, list(map(int, vec_gamma))).compute(), Z_u_a), SingleMSM(self.vec_R_U, list(map(int, vec_gamma_inv))).compute())
    msm_accumulator.accumulate_check(point_lhs, list(map(affine_to_projective, vec_U)), list(map(int, vec_x_times_s)))

    return True, ''

def test_same_msm():
  transcript_prover = CurdleproofsTranscript()
  n = 128

  crs_G_vec = [normalize(get_random_point()) for _ in range(0, n)]

  vec_T = [normalize(get_random_point()) for _ in range(0, n)]
  vec_U = [normalize(get_random_point()) for _ in range(0, n)]
  vec_x = [Fr(random.randint(1, Fr.field_modulus - 1)) for _ in range(0, n)]

  A = SingleMSM(list(map(affine_to_projective, crs_G_vec)), list(map(int, vec_x))).compute()
  Z_t = SingleMSM(list(map(affine_to_projective, vec_T)), list(map(int, vec_x))).compute()
  Z_u = SingleMSM(list(map(affine_to_projective, vec_U)), list(map(int, vec_x))).compute()

  proof: SameMSMProof = SameMSMProof.new(
    crs_G_vec=crs_G_vec,
    A=A,
    Z_t=Z_t,
    Z_u=Z_u,
    vec_T=vec_T,
    vec_U=vec_U,
    vec_x=vec_x,
    transcript=transcript_prover
  )

  print("Proof", proof)

  transcript_verifier = CurdleproofsTranscript()
  msm_accumulator = MSMAccumulatorInefficient()

  (result, err) = proof.verify(
    crs_G_vec=list(map(affine_to_projective, crs_G_vec)),
    A=A,
    Z_t=Z_t,
    Z_u=Z_u,
    vec_T=vec_T,
    vec_U=vec_U,
    transcript=transcript_verifier,
    msm_accumulator=msm_accumulator
  )

  msm_verify = msm_accumulator.verify()
  print("Result", result)
  print("MSM verify", msm_verify)
  print("Error", err)

test_same_msm()