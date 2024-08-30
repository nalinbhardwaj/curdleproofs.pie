from functools import reduce
import operator
import random
import pytest
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.grand_prod import GrandProductProof
from curdleproofs.opening import TrackerOpeningProof
from curdleproofs.util import get_random_point, get_permutation
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from typing import List
from curdleproofs.util import (
    generate_blinders,
    inner_product,
    field_to_bytes,
    point_projective_to_bytes,
    random_scalar,
    scalar_pow,
    G1,
    Z1,
    CURVE_ORDER,
    BLSPubkey,
)
from curdleproofs.msm_accumulator import MSMAccumulator, compute_MSM
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
from py_arkworks_bls12381 import G1Point, Scalar


def test_py_arkworks_bls12381_api():
    print(dir(G1Point))
    assert dir(G1Point) == [
        '__add__',
        '__class__',
        '__delattr__',
        '__dir__',
        '__doc__',
        '__eq__',
        '__format__',
        '__ge__',
        '__getattribute__',
        '__gt__',
        '__hash__',
        '__init__',
        '__init_subclass__',
        '__le__',
        '__lt__',
        '__module__',
        '__mul__',
        '__ne__',
        '__neg__',
        '__new__',
        '__radd__',
        '__reduce__',
        '__reduce_ex__',
        '__repr__',
        '__rmul__',
        '__rsub__',
        '__setattr__',
        '__sizeof__',
        '__str__',
        '__sub__',
        '__subclasshook__',
        'from_compressed_bytes',
        'from_compressed_bytes_unchecked',
        'identity',
        'multiexp_unchecked',
        'to_compressed_bytes'
    ]

    print(dir(Scalar))
    assert dir(Scalar) == [
        '__add__',
        '__class__',
        '__delattr__',
        '__dir__',
        '__doc__',
        '__eq__',
        '__format__',
        '__ge__',
        '__getattribute__',
        '__gt__',
        '__hash__',
        '__init__',
        '__init_subclass__',
        '__int__',
        '__le__',
        '__lt__',
        '__module__',
        '__mul__',
        '__ne__',
        '__neg__',
        '__new__',
        '__radd__',
        '__reduce__',
        '__reduce_ex__',
        '__repr__',
        '__rmul__',
        '__rsub__',
        '__rtruediv__',
        '__setattr__',
        '__sizeof__',
        '__str__',
        '__sub__',
        '__subclasshook__',
        '__truediv__',
        'from_le_bytes',
        'inverse',
        'is_zero',
        'pow',
        'square',
        'to_le_bytes'
    ]


# Copied from https://pypi.org/project/py-arkworks-bls12381/
def test_py_arkworks_bls12381_g1points():
    # G1Point and G2Point have the same methods implemented on them
    # For brevity, I will only show one method using G1Point and G2Point
    # The rest of the code will just use G1Point

    # Point initialization -- This will be initialized to the g1 generator
    g1_generator = G1Point()

    # Identity element
    identity = G1Point.identity()

    # Equality -- We override eq and neq operators
    assert g1_generator == g1_generator
    assert g1_generator != identity

    # Printing an element -- We override __str__ so when we print
    # an element it prints in hex
    print("identity: ", identity)
    print("g1 generator: ", g1_generator)

    # Point Addition/subtraction/Negation -- We override the add/sub/neg operators
    gen = G1Point()
    double_gen = gen + gen
    assert double_gen - gen == gen
    neg_gen = -gen
    assert neg_gen + gen == identity

    # Scalar multiplication
    #
    scalar = Scalar(4)
    four_gen = gen * scalar
    assert four_gen == gen + gen + gen + gen

    # Serialisation
    #
    # serialising to/from a g1 point
    # We don't expose the uncompressed form
    # because it seems like its not needed
    compressed_bytes = gen.to_compressed_bytes()
    deserialised_point = G1Point.from_compressed_bytes(compressed_bytes)
    # If the bytes being received are trusted, we can avoid
    # doing subgroup checks
    deserialised_point_unchecked = G1Point.from_compressed_bytes_unchecked(compressed_bytes)
    assert deserialised_point == deserialised_point_unchecked
    assert deserialised_point == gen

    # Serialization
    assert str(gen) == "97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb"
    assert point_projective_to_bytes(gen) == bytes.fromhex("97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb")

    # Indexing
    point_a = G1 * Scalar(4)
    point_a_copy = G1 * Scalar(4)
    point_map = {}
    # G1Point should not be hashable
    with pytest.raises(TypeError):
        point_map[point_a] = True
    # Able to index by serialized form (inefficient)
    point_map[point_projective_to_bytes(point_a)] = True
    assert point_map[point_projective_to_bytes(point_a_copy)]


def test_py_arkworks_bls12381_scalar():
    scalar = Scalar(4)
    assert field_to_bytes(scalar) == bytes.fromhex("0400000000000000000000000000000000000000000000000000000000000000")

    assert CURVE_ORDER == 52435875175126190479447740508185965837690552500527637822603658699938581184513

    # Scalar should handle Fr::MAX
    assert int(Scalar(CURVE_ORDER - 1)) == CURVE_ORDER - 1
    # Scalar should overflow if given a value >= CURVE_ORDER
    assert int(Scalar(CURVE_ORDER)) == 0
    # Sanity check that the overflowed value is correct
    assert int(Scalar(2**256)) == 2**256 % CURVE_ORDER
    # It even accepts integers with more than 256 bits
    assert int(Scalar(2**257)) == 2**257 % CURVE_ORDER

    # Deserialize from big integers
    Scalar.from_le_bytes((CURVE_ORDER - 1).to_bytes(32, 'little'))
    with pytest.raises(ValueError):
        # Errors with `ValueError: Err From Rust: serialised data seems to be invalid`
        Scalar.from_le_bytes((CURVE_ORDER).to_bytes(32, 'little'))


def test_scalar_pow():
    helper_test_scalar_pow(1, 1)
    helper_test_scalar_pow(4, 1)
    helper_test_scalar_pow(4, 2)
    helper_test_scalar_pow(4, 3)
    helper_test_scalar_pow(4, 4)
    helper_test_scalar_pow(100, 2)
    helper_test_scalar_pow(42, 6)
    scalar_pow(random_scalar(), 128)


def helper_test_scalar_pow(base: int, exponent: int):
    res_scalar = scalar_pow(Scalar(base), exponent)
    res = int.from_bytes(res_scalar.to_le_bytes(), byteorder='little')
    assert res == base ** exponent


def test_utils_point_projective_to_bytes():
    scalar = Scalar(99)
    point = G1 * scalar
    point_projective_to_bytes(point) == bytes.fromhex("aa10e1055b14a89cc3261699524998732fddc4f30c76c1057eb83732a01416643eb015a932e4080c86f42e485973d240")


def test_utils_get_random_point():
    point = get_random_point()
    assert point + G1 - point == G1


def test_ipa():
    transcript = CurdleproofsTranscript(b"curdleproofs")

    n = 128

    crs_G_vec = [get_random_point() for _ in range(0, n)]

    vec_u = generate_blinders(n)
    crs_G_prime_vec = [G_i * u_i for (G_i, u_i) in zip(crs_G_vec, vec_u)]
    crs_H = get_random_point()

    vec_b = [random_scalar() for _ in range(0, n)]
    vec_c = [random_scalar() for _ in range(0, n)]

    z = inner_product(vec_b, vec_c)

    B = compute_MSM(crs_G_vec, vec_b)
    C = compute_MSM(crs_G_prime_vec, vec_c)

    proof = IPA.new(
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

    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator = MSMAccumulator()

    proof.verify(
        crs_G_vec=crs_G_vec,
        crs_H=crs_H,
        C=B,
        D=C,
        inner_prod=z,
        vec_u=vec_u,
        transcript=transcript_verifier,
        msm_accumulator=msm_accumulator,
    )
    msm_accumulator.verify()

    transcript_wrong = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator_wrong = MSMAccumulator()

    with pytest.raises(AssertionError):
        proof.verify(
            crs_G_vec=crs_G_vec,
            crs_H=crs_H,
            C=B,
            D=C,
            inner_prod=z + Scalar(1),
            vec_u=vec_u,
            transcript=transcript_wrong,
            msm_accumulator=msm_accumulator_wrong,
        )
        msm_accumulator_wrong.verify()


def test_gprod():
    transcript_prover = CurdleproofsTranscript(b"curdleproofs")

    n = 128
    n_blinders = 4
    ell = n - n_blinders

    crs_G_vec = [get_random_point() for _ in range(ell)]
    crs_H_vec = [get_random_point() for _ in range(n_blinders)]
    crs_U = get_random_point()
    crs_G_sum = reduce(lambda a, b: a + b, crs_G_vec, Z1)
    crs_H_sum = reduce(lambda a, b: a + b, crs_H_vec, Z1)

    vec_b = [random_scalar() for _ in range(ell)]
    vec_b_blinders = generate_blinders(n_blinders)

    gprod_result = reduce(operator.mul, vec_b, Scalar(1))

    B = compute_MSM(crs_G_vec, vec_b) + compute_MSM(crs_H_vec, vec_b_blinders)

    gprod_proof = GrandProductProof.new(
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

    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator = MSMAccumulator()

    gprod_proof.verify(
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

    msm_accumulator.verify()

    # Wrong test
    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator = MSMAccumulator()
    with pytest.raises(AssertionError):
        gprod_proof.verify(
            crs_G_vec=crs_G_vec,
            crs_H_vec=crs_H_vec,
            crs_U=crs_U,
            crs_G_sum=crs_G_sum,
            crs_H_sum=crs_H_sum,
            B=B,
            gprod_result=gprod_result + Scalar(1),
            n_blinders=n_blinders,
            transcript=transcript_verifier,
            msm_accumulator=msm_accumulator,
        )

        msm_accumulator.verify()

    # Wrong test
    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator = MSMAccumulator()
    with pytest.raises(AssertionError):
        gprod_proof.verify(
            crs_G_vec=crs_G_vec,
            crs_H_vec=crs_H_vec,
            crs_U=crs_U,
            crs_G_sum=crs_G_sum,
            crs_H_sum=crs_H_sum,
            B=B * Scalar(3),
            gprod_result=gprod_result,
            n_blinders=n_blinders,
            transcript=transcript_verifier,
            msm_accumulator=msm_accumulator,
        )

        msm_accumulator.verify()


def test_same_permutation_proof():
    transcript_prover = CurdleproofsTranscript(b"curdleproofs")

    n = 128
    n_blinders = 4
    ell = n - n_blinders

    crs_G_vec = [get_random_point() for _ in range(0, ell)]
    crs_H_vec = [get_random_point() for _ in range(0, n_blinders)]

    crs_U = get_random_point()
    crs_G_sum = reduce(lambda a, b: a + b, crs_G_vec, Z1)
    crs_H_sum = reduce(lambda a, b: a + b, crs_H_vec, Z1)

    vec_a_blinders = generate_blinders(n_blinders)
    vec_m_blinders = generate_blinders(n_blinders)

    permutation = list(range(0, ell))
    random.shuffle(permutation)

    vec_a = [random_scalar() for _ in range(0, ell)]
    vec_a_permuted = get_permutation(vec_a, permutation)

    A = compute_MSM(crs_G_vec, vec_a_permuted) + compute_MSM(crs_H_vec, vec_a_blinders)
    M = compute_MSM(crs_G_vec, map(Scalar, permutation)) + compute_MSM(crs_H_vec, vec_m_blinders)

    same_perm_proof = SamePermutationProof.new(
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

    transcript_verifier = CurdleproofsTranscript(b"curdleproofs")
    msm_accumulator = MSMAccumulator()

    same_perm_proof.verify(
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

    msm_accumulator.verify()


def test_same_msm():
    transcript_prover = CurdleproofsTranscript(b"curdleproofs")
    n = 128

    crs_G_vec = [get_random_point() for _ in range(0, n)]

    vec_T = [get_random_point() for _ in range(0, n)]
    vec_U = [get_random_point() for _ in range(0, n)]
    vec_x = [random_scalar() for _ in range(0, n)]

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

    proof.verify(
        crs_G_vec=crs_G_vec,
        A=A,
        Z_t=Z_t,
        Z_u=Z_u,
        vec_T=vec_T,
        vec_U=vec_U,
        transcript=transcript_verifier,
        msm_accumulator=msm_accumulator,
    )

    msm_accumulator.verify()


def test_same_scalar_arg():
    transcript_prover = CurdleproofsTranscript(b"curdleproofs")

    crs_G_t = get_random_point()
    crs_G_u = get_random_point()
    crs_H = get_random_point()

    R = get_random_point()
    S = get_random_point()

    k = random_scalar()
    r_t = random_scalar()
    r_u = random_scalar()

    cm_T = GroupCommitment.new(crs_G_t, crs_H, R * k, r_t)
    cm_U = GroupCommitment.new(crs_G_u, crs_H, S * k, r_u)

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
    proof.verify(
        crs_G_t=crs_G_t,
        crs_G_u=crs_G_u,
        crs_H=crs_H,
        R=R,
        S=S,
        cm_T=cm_T,
        cm_U=cm_U,
        transcript=transcript_verifier,
    )


def test_group_commit():
    crs_G = get_random_point()
    crs_H = get_random_point()

    A = get_random_point()
    B = get_random_point()

    r_a = random_scalar()
    r_b = random_scalar()

    cm_a = GroupCommitment.new(crs_G, crs_H, A, r_a)
    cm_b = GroupCommitment.new(crs_G, crs_H, B, r_b)
    cm_a_b = GroupCommitment.new(crs_G, crs_H, A + B, r_a + r_b)

    assert cm_a + cm_b == cm_a_b


def test_shuffle_argument():
    N = 64
    ell = N - N_BLINDERS

    crs = CurdleproofsCrs.new(ell, N_BLINDERS)

    permutation = list(range(ell))
    random.shuffle(permutation)
    k = random_scalar()

    vec_R = [get_random_point() for _ in range(ell)]
    vec_S = [get_random_point() for _ in range(ell)]

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
    shuffle_proof.verify(crs, vec_R, vec_S, vec_T, vec_U, M)


def test_bad_shuffle_argument():
    N = 128
    ell = N - N_BLINDERS

    crs = CurdleproofsCrs.new(ell, N_BLINDERS)

    permutation = list(range(ell))
    random.shuffle(permutation)
    k = random_scalar()

    vec_R = [get_random_point() for _ in range(ell)]
    vec_S = [get_random_point() for _ in range(ell)]

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

    with pytest.raises(AssertionError):
        shuffle_proof.verify(crs, vec_S, vec_R, vec_T, vec_U, M)

    another_permutation = list(range(ell))
    random.shuffle(another_permutation)

    with pytest.raises(AssertionError):
        shuffle_proof.verify(
            crs,
            vec_R,
            vec_S,
            get_permutation(vec_T, another_permutation),
            get_permutation(vec_U, another_permutation),
            M,
        )

        shuffle_proof.verify(
            crs, vec_R, vec_S, vec_T, vec_U, M * k
        )

    another_k = random_scalar()
    another_vec_T = [T * another_k for T in vec_T]
    another_vec_U = [U * another_k for U in vec_U]

    with pytest.raises(AssertionError):
        shuffle_proof.verify(
            crs, vec_R, vec_S, another_vec_T, another_vec_U, M
        )


def test_serde():
    N = 64
    ell = N - N_BLINDERS

    crs = CurdleproofsCrs.new(ell, N_BLINDERS)

    permutation = list(range(ell))
    random.shuffle(permutation)
    k = random_scalar()

    vec_R = [get_random_point() for _ in range(ell)]
    vec_S = [get_random_point() for _ in range(ell)]

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
    deser_shuffle_proof.verify(
        deser_crs,
        deser_verifier_input.vec_R,
        deser_verifier_input.vec_S,
        deser_verifier_input.vec_T,
        deser_verifier_input.vec_U,
        deser_verifier_input.M,
    )


def test_tracker_opening_proof():
    G = G1
    k = generate_blinders(1)[0]
    r = generate_blinders(1)[0]

    k_G = G * k
    r_G = G * r
    k_r_G = r_G * k

    transcript_prover = CurdleproofsTranscript(b"whisk_opening_proof")
    opening_proof = TrackerOpeningProof.new(
        k_r_G=k_r_G, r_G=r_G, k_G=k_G, k=k, transcript=transcript_prover
    )

    json_str_proof = opening_proof.to_json()
    print("json_str_proof", json_str_proof)

    deser_proof = TrackerOpeningProof.from_json(json_str_proof)

    transcript_verifier = CurdleproofsTranscript(b"whisk_opening_proof")
    deser_proof.verify(transcript_verifier, k_r_G, r_G, k_G)


def test_whisk_interface_tracker_opening_proof():
    k = generate_random_k()
    k_commitment = get_k_commitment(k)
    tracker = generate_tracker(k)

    tracker_proof = GenerateWhiskTrackerProof(tracker, k)

    IsValidWhiskOpeningProof(tracker, k_commitment, tracker_proof)


def test_whisk_interface_shuffle_proof():
    N = 64
    ell = N - N_BLINDERS
    crs = generate_random_crs(ell)
    pre_trackers = generate_random_trackers(ell)
    post_trackers, shuffle_proof = GenerateWhiskShuffleProof(crs, pre_trackers)
    IsValidWhiskShuffleProof(crs, pre_trackers, post_trackers, shuffle_proof)


def generate_random_k() -> Scalar:
    return generate_blinders(1)[0]


def get_k_commitment(k: Scalar) -> BLSPubkey:
    return BLSPubkey(point_projective_to_bytes(G1 * k))


def generate_tracker(k: Scalar) -> WhiskTracker:
    r = generate_blinders(1)[0]
    r_G = G1 * r
    k_r_G = r_G * k
    return WhiskTracker(BLSPubkey(point_projective_to_bytes(r_G)), BLSPubkey(point_projective_to_bytes(k_r_G)))


def generate_random_crs(ell: int) -> CurdleproofsCrs:
    return CurdleproofsCrs.new(ell, N_BLINDERS)


def generate_random_trackers(n: int) -> List[WhiskTracker]:
    return [generate_tracker(generate_random_k()) for _ in range(n)]
