from functools import reduce
from math import log2
import operator
import random
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.grand_prod import GrandProductProof
from curdleproofs.opening import TrackerOpeningProof
from curdleproofs.util import (
    affine_to_projective,
    point_affine_to_bytes,
    point_projective_to_bytes,
    points_affine_to_bytes,
    points_projective_to_bytes,
    get_random_point,
    get_permutation,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List, Optional, Tuple, Type, TypeVar
from curdleproofs.util import (
    PointAffine,
    PointProjective,
    Fr,
    field_to_bytes,
    invert,
    generate_blinders,
    inner_product,
    get_verification_scalars_bitstring,
)
from curdleproofs.msm_accumulator import MSMAccumulator, compute_MSM
from py_ecc.optimized_bls12_381.optimized_curve import (
    curve_order,
    G1,
    multiply,
    normalize,
    add,
    neg,
    Z1,
)
from curdleproofs.ipa import IPA
from curdleproofs.same_perm import SamePermutationProof
from curdleproofs.same_msm import SameMSMProof
from curdleproofs.commitment import GroupCommitment
from curdleproofs.same_scalar import SameScalarProof
from curdleproofs.curdleproofs import (
    N_BLINDERS,
    CurdleProofsProof,
    VerifierInput,
    shuffle_permute_and_commit_input,
)
from curdleproofs.whisk_interface import (
    WhiskTracker,
    GenerateWhiskTrackerProof,
    GenerateWhiskShuffleProof,
    IsValidWhiskOpeningProof,
    IsValidWhiskShuffleProof,
)


def test_ipa():
    transcript = CurdleproofsTranscript(b"curdleproofs")

    n = 128

    crs_G_vec = [get_random_point() for _ in range(0, n)]

    vec_u = generate_blinders(n)
    crs_G_prime_vec = [multiply(G_i, int(u_i)) for (G_i, u_i) in zip(crs_G_vec, vec_u)]
    crs_H = get_random_point()

    vec_b = [Fr(random.randint(1, Fr.field_modulus)) for _ in range(0, n)]
    vec_c = [Fr(random.randint(1, Fr.field_modulus)) for _ in range(0, n)]

    z = inner_product(vec_b, vec_c)
    print("prod = ", vec_b, vec_c, z)

    B = compute_MSM(crs_G_vec, vec_b)
    C = compute_MSM(crs_G_prime_vec, vec_c)

    (proof, err) = IPA.new(
        crs_G_vec=crs_G_vec,
        crs_G_prime_vec=crs_G_prime_vec,
        crs_H=crs_H,
        C=B,
        D=C,
        z=z,
        vec_c=vec_b,
        vec_d=vec_c,
        transcript=transcript,
    )

    print(
        "proof: ",
        proof,
        proof.vec_L_C,
        proof.vec_L_D,
        proof.vec_R_C,
        proof.vec_R_D,
        "crs len",
        len(crs_G_vec),
    )
    print("err: ", err)

    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator = MSMAccumulator()

    (result, err) = proof.verify(
        crs_G_vec=crs_G_vec,
        crs_H=crs_H,
        C=B,
        D=C,
        inner_prod=z,
        vec_u=vec_u,
        transcript=transcript_verifier,
        msm_accumulator=msm_accumulator,
    )
    msm_verify = msm_accumulator.verify()

    print("result: ", result)
    print("msm_verify: ", msm_verify)
    print("err: ", err)
    assert result and msm_verify

    transcript_wrong = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator_wrong = MSMAccumulator()
    (result_wrong, err_wrong) = proof.verify(
        crs_G_vec=crs_G_vec,
        crs_H=crs_H,
        C=B,
        D=C,
        inner_prod=z + Fr.one(),
        vec_u=vec_u,
        transcript=transcript_wrong,
        msm_accumulator=msm_accumulator_wrong,
    )
    msm_wrong_verify = msm_accumulator_wrong.verify()
    print("result_wrong: ", result_wrong)
    print("msm_wrong_verify: ", msm_wrong_verify)
    print("err_wrong: ", err_wrong)
    assert not (result_wrong and msm_wrong_verify)


def test_gprod():
    transcript_prover = CurdleproofsTranscript(b"curdleproofs")

    n = 128
    n_blinders = 4
    ell = n - n_blinders

    crs_G_vec = [get_random_point() for _ in range(ell)]
    crs_H_vec = [get_random_point() for _ in range(n_blinders)]
    crs_U = get_random_point()
    crs_G_sum = reduce(add, crs_G_vec, Z1)
    crs_H_sum = reduce(add, crs_H_vec, Z1)

    vec_b = [Fr(random.randint(1, Fr.field_modulus - 1)) for _ in range(ell)]
    vec_b_blinders = generate_blinders(n_blinders)

    gprod_result = reduce(operator.mul, vec_b, Fr.one())

    B = add(compute_MSM(crs_G_vec, vec_b), compute_MSM(crs_H_vec, vec_b_blinders))

    (gprod_proof, err) = GrandProductProof.new(
        crs_G_vec=crs_G_vec,
        crs_H_vec=crs_H_vec,
        crs_U=crs_U,
        B=B,
        gprod_result=gprod_result,
        vec_b=vec_b,
        vec_b_blinders=vec_b_blinders,
        transcript=transcript_prover,
    )

    print("Prover result: ", gprod_proof)
    print("Prover error:", err)

    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator = MSMAccumulator()

    (result, err) = gprod_proof.verify(
        crs_G_vec=crs_G_vec,
        crs_H_vec=crs_H_vec,
        crs_U=crs_U,
        crs_G_sum=crs_G_sum,
        crs_H_sum=crs_H_sum,
        B=B,
        gprod_result=gprod_result,
        n_blinders=n_blinders,
        transcript=transcript_verifier,
        msm_accumulator=msm_accumulator,
    )

    msm_verify = msm_accumulator.verify()

    print("Result: ", result)
    print("MSM verify: ", msm_verify)
    print("Error: ", err)
    assert result and msm_verify

    # Wrong test
    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator = MSMAccumulator()
    (result, err) = gprod_proof.verify(
        crs_G_vec=crs_G_vec,
        crs_H_vec=crs_H_vec,
        crs_U=crs_U,
        crs_G_sum=crs_G_sum,
        crs_H_sum=crs_H_sum,
        B=B,
        gprod_result=gprod_result + Fr.one(),
        n_blinders=n_blinders,
        transcript=transcript_verifier,
        msm_accumulator=msm_accumulator,
    )

    msm_verify = msm_accumulator.verify()

    print("Wrong Result: ", result)
    print("Wrong MSM verify: ", msm_verify)
    print("Wrong Error: ", err)
    assert not (result and msm_verify)

    # Wrong test
    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator = MSMAccumulator()
    (result, err) = gprod_proof.verify(
        crs_G_vec=crs_G_vec,
        crs_H_vec=crs_H_vec,
        crs_U=crs_U,
        crs_G_sum=crs_G_sum,
        crs_H_sum=crs_H_sum,
        B=multiply(B, 3),
        gprod_result=gprod_result,
        n_blinders=n_blinders,
        transcript=transcript_verifier,
        msm_accumulator=msm_accumulator,
    )

    msm_verify = msm_accumulator.verify()

    print("Wrong Result: ", result)
    print("Wrong MSM verify: ", msm_verify)
    print("Wrong Error: ", err)
    assert not (result and msm_verify)


def test_same_permutation_proof():
    transcript_prover = CurdleproofsTranscript(b"curdleproofs")

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

    A = add(
        compute_MSM(crs_G_vec, vec_a_permuted), compute_MSM(crs_H_vec, vec_a_blinders)
    )
    M = add(compute_MSM(crs_G_vec, permutation), compute_MSM(crs_H_vec, vec_m_blinders))

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
        transcript=transcript_prover,
    )

    print("Proof: ", same_perm_proof)
    print("Error: ", err)

    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator = MSMAccumulator()

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
        msm_accumulator=msm_accumulator,
    )

    msm_verify = msm_accumulator.verify()

    print("Verify: ", verify)
    print("Error: ", err)
    print("MSM verify: ", msm_verify)
    assert verify and msm_verify


def test_same_msm():
    transcript_prover = CurdleproofsTranscript(b"curdleproofs")
    n = 128

    crs_G_vec = [get_random_point() for _ in range(0, n)]

    vec_T = [get_random_point() for _ in range(0, n)]
    vec_U = [get_random_point() for _ in range(0, n)]
    vec_x = [Fr(random.randint(1, Fr.field_modulus)) for _ in range(0, n)]

    A = compute_MSM(crs_G_vec, vec_x)
    Z_t = compute_MSM(vec_T, vec_x)
    Z_u = compute_MSM(vec_U, vec_x)

    proof: SameMSMProof = SameMSMProof.new(
        crs_G_vec=crs_G_vec,
        A=A,
        Z_t=Z_t,
        Z_u=Z_u,
        vec_T=vec_T,
        vec_U=vec_U,
        vec_x=vec_x,
        transcript=transcript_prover,
    )

    print("Proof", proof)

    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator = MSMAccumulator()

    (result, err) = proof.verify(
        crs_G_vec=crs_G_vec,
        A=A,
        Z_t=Z_t,
        Z_u=Z_u,
        vec_T=vec_T,
        vec_U=vec_U,
        transcript=transcript_verifier,
        msm_accumulator=msm_accumulator,
    )

    msm_verify = msm_accumulator.verify()
    print("Result", result)
    print("MSM verify", msm_verify)
    print("Error", err)
    assert result and msm_verify


def test_same_scalar_arg():
    transcript_prover = CurdleproofsTranscript(b"curdleproofs")

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
        transcript=transcript_prover,
    )

    print("proof", proof)

    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    (res, err) = proof.verify(
        crs_G_t=crs_G_t,
        crs_G_u=crs_G_u,
        crs_H=crs_H,
        R=R,
        S=S,
        cm_T=cm_T,
        cm_U=cm_U,
        transcript=transcript_verifier,
    )
    print("res", res)
    print("err", err)
    assert res


def test_group_commit():
    crs_G = get_random_point()
    crs_H = get_random_point()

    A = get_random_point()
    B = get_random_point()

    r_a = Fr(random.randint(1, Fr.field_modulus))
    r_b = Fr(random.randint(1, Fr.field_modulus))

    cm_a = GroupCommitment.new(crs_G, crs_H, A, r_a)
    cm_b = GroupCommitment.new(crs_G, crs_H, B, r_b)
    cm_a_b = GroupCommitment.new(crs_G, crs_H, add(A, B), r_a + r_b)

    assert cm_a + cm_b == cm_a_b


def test_shuffle_argument():
    N = 64
    ell = N - N_BLINDERS

    crs = CurdleproofsCrs.new(ell, N_BLINDERS)

    permutation = list(range(ell))
    random.shuffle(permutation)
    k = Fr(random.randint(1, Fr.field_modulus))

    vec_R = [normalize(get_random_point()) for _ in range(ell)]
    vec_S = [normalize(get_random_point()) for _ in range(ell)]

    vec_T, vec_U, M, vec_m_blinders = shuffle_permute_and_commit_input(
        crs, vec_R, vec_S, permutation, k
    )

    shuffle_proof: CurdleProofsProof = CurdleProofsProof.new(
        crs=crs,
        vec_R=vec_R,
        vec_S=vec_S,
        vec_T=vec_T,
        vec_U=vec_U,
        M=M,
        permutation=permutation,
        k=k,
        vec_m_blinders=vec_m_blinders,
    )

    print("shuffle proof", shuffle_proof)

    # for i in range(50):
    #   print("iter ", i)
    verify, err = shuffle_proof.verify(crs, vec_R, vec_S, vec_T, vec_U, M)
    print("verify", verify)
    print("err", err)
    assert verify


def test_bad_shuffle_argument():
    N = 128
    ell = N - N_BLINDERS

    crs = CurdleproofsCrs.new(ell, N_BLINDERS)

    permutation = list(range(ell))
    random.shuffle(permutation)
    k = Fr(random.randint(1, Fr.field_modulus))

    vec_R = [normalize(get_random_point()) for _ in range(ell)]
    vec_S = [normalize(get_random_point()) for _ in range(ell)]

    vec_T, vec_U, M, vec_m_blinders = shuffle_permute_and_commit_input(
        crs, vec_R, vec_S, permutation, k
    )

    shuffle_proof: CurdleProofsProof = CurdleProofsProof.new(
        crs=crs,
        vec_R=vec_R,
        vec_S=vec_S,
        vec_T=vec_T,
        vec_U=vec_U,
        M=M,
        permutation=permutation,
        k=k,
        vec_m_blinders=vec_m_blinders,
    )

    print("shuffle proof", shuffle_proof)

    verify, err = shuffle_proof.verify(crs, vec_S, vec_R, vec_T, vec_U, M)
    print("false verify", verify)
    print("err", err)
    assert not verify

    another_permutation = list(range(ell))
    random.shuffle(another_permutation)

    verify, err = shuffle_proof.verify(
        crs,
        vec_R,
        vec_S,
        get_permutation(vec_T, another_permutation),
        get_permutation(vec_U, another_permutation),
        M,
    )
    print("false verify also", verify)
    print("err", err)
    assert not verify

    verify, err = shuffle_proof.verify(
        crs, vec_R, vec_S, vec_T, vec_U, multiply(M, int(k))
    )
    print("false verify also also", verify)
    print("err", err)
    assert not verify

    another_k = Fr(random.randint(1, Fr.field_modulus))
    another_vec_T = [
        normalize(multiply(affine_to_projective(T), int(another_k))) for T in vec_T
    ]
    another_vec_U = [
        normalize(multiply(affine_to_projective(U), int(another_k))) for U in vec_U
    ]

    verify, err = shuffle_proof.verify(
        crs, vec_R, vec_S, another_vec_T, another_vec_U, M
    )
    print("false verify also also also", verify)
    print("err", err)
    assert not verify


def test_serde():
    N = 64
    ell = N - N_BLINDERS

    crs = CurdleproofsCrs.new(ell, N_BLINDERS)

    permutation = list(range(ell))
    random.shuffle(permutation)
    k = Fr(random.randint(1, Fr.field_modulus))

    vec_R = [normalize(get_random_point()) for _ in range(ell)]
    vec_S = [normalize(get_random_point()) for _ in range(ell)]

    vec_T, vec_U, M, vec_m_blinders = shuffle_permute_and_commit_input(
        crs, vec_R, vec_S, permutation, k
    )

    shuffle_proof: CurdleProofsProof = CurdleProofsProof.new(
        crs=crs,
        vec_R=vec_R,
        vec_S=vec_S,
        vec_T=vec_T,
        vec_U=vec_U,
        M=M,
        permutation=permutation,
        k=k,
        vec_m_blinders=vec_m_blinders,
    )

    print("shuffle proof", shuffle_proof)

    json_str_proof = shuffle_proof.to_json()
    print("json_str_proof", json_str_proof)

    json_str_crs = crs.to_json()
    print("json_str_crs", json_str_crs)

    verifier_input = VerifierInput(
        vec_R=vec_R,
        vec_S=vec_S,
        vec_T=vec_T,
        vec_U=vec_U,
        M=M,
    )

    json_str_verifier_input = verifier_input.to_json()

    deser_shuffle_proof = CurdleProofsProof.from_json(json_str_proof)
    deser_crs = CurdleproofsCrs.from_json(json_str_crs)
    deser_verifier_input = VerifierInput.from_json(json_str_verifier_input)

    # for i in range(50):
    #   print("iter ", i)
    verify, err = deser_shuffle_proof.verify(
        deser_crs,
        deser_verifier_input.vec_R,
        deser_verifier_input.vec_S,
        deser_verifier_input.vec_T,
        deser_verifier_input.vec_U,
        deser_verifier_input.M,
    )
    print("verify", verify)
    print("err", err)
    assert verify


def test_tracker_opening_proof():
    G = G1
    k = generate_blinders(1)[0]
    r = generate_blinders(1)[0]

    k_G = multiply(G, int(k))
    r_G = multiply(G, int(r))
    k_r_G = multiply(r_G, int(k))

    transcript_prover = CurdleproofsTranscript(b"whisk_opening_proof")
    opening_proof = TrackerOpeningProof.new(
        k_r_G=k_r_G, r_G=r_G, k_G=k_G, G=G, k=k, transcript=transcript_prover
    )

    json_str_proof = opening_proof.to_json()
    print("json_str_proof", json_str_proof)

    deser_proof = TrackerOpeningProof.from_json(json_str_proof)

    transcript_verifier = CurdleproofsTranscript(b"whisk_opening_proof")
    assert deser_proof.verify(transcript_verifier, k_r_G, r_G, k_G)


def test_whisk_interface_tracker_opening_proof():
    k = generate_random_k()
    tracker = generate_tracker(k)

    k_commitment = get_k_commitment(k)
    tracker_proof = GenerateWhiskTrackerProof(tracker, k_commitment, k)

    assert IsValidWhiskOpeningProof(tracker, k_commitment, tracker_proof)


def test_whisk_interface_shuffle_proof():
    N = 64
    ell = N - N_BLINDERS
    crs = generate_random_crs(ell)
    pre_trackers = generate_random_trackers(ell)
    post_trackers, m, shuffle_proof = GenerateWhiskShuffleProof(crs, pre_trackers)
    assert IsValidWhiskShuffleProof(crs, pre_trackers, post_trackers, m, shuffle_proof)


def generate_random_k() -> Fr:
    return generate_blinders(1)[0]


def get_k_commitment(k: Fr) -> PointAffine:
    return normalize(multiply(G1, int(k)))


def generate_tracker(k: Fr) -> WhiskTracker:
    r = generate_blinders(1)[0]
    r_G = multiply(G1, int(r))
    k_r_G = multiply(r_G, int(k))
    return WhiskTracker(normalize(r_G), normalize(k_r_G))


def generate_random_crs(ell: int) -> CurdleproofsCrs:
    return CurdleproofsCrs.new(ell, N_BLINDERS)


def generate_random_trackers(n: int) -> List[WhiskTracker]:
    return [generate_tracker(generate_random_k()) for _ in range(n)]
