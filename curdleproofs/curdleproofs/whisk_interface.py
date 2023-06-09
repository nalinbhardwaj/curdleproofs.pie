import random
from typing import Container, Sequence, Tuple
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.curdleproofs import (
    N_BLINDERS,
    CurdleProofsProof,
    shuffle_permute_and_commit_input,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from curdleproofs.opening import TrackerOpeningProof
from curdleproofs.util import (
    PointAffine as BLSG1Point,
    affine_to_projective,
    Fr,
)
from py_ecc.optimized_bls12_381.optimized_curve import G1, normalize


class WhiskTracker:
    r_G: BLSG1Point  # r * G
    k_r_G: BLSG1Point  # k * r * G

    def __init__(self, r_G: BLSG1Point, k_r_G: BLSG1Point):
        self.r_G = r_G
        self.k_r_G = k_r_G


SerializedCurdleProofsProof = bytes


def IsValidWhiskShuffleProof(
    crs: CurdleproofsCrs,
    pre_shuffle_trackers: Sequence[WhiskTracker],
    post_shuffle_trackers: Sequence[WhiskTracker],
    shuffle_proof: SerializedCurdleProofsProof,
) -> Tuple[bool, str]:
    """
    Verify `post_shuffle_trackers` is a permutation of `pre_shuffle_trackers`.
    """
    vec_R = [tracker.r_G for tracker in pre_shuffle_trackers]
    vec_S = [tracker.k_r_G for tracker in pre_shuffle_trackers]

    vec_T = [tracker.r_G for tracker in post_shuffle_trackers]
    vec_U = [tracker.k_r_G for tracker in post_shuffle_trackers]

    shuffle_proof_instance = CurdleProofsProof.from_json(shuffle_proof.decode())

    return shuffle_proof_instance.verify(crs, vec_R, vec_S, vec_T, vec_U)


def GenerateWhiskShuffleProof(
    crs: CurdleproofsCrs, pre_shuffle_trackers: Sequence[WhiskTracker]
) -> Tuple[Sequence[WhiskTracker], SerializedCurdleProofsProof]:
    permutation = list(range(len(crs.vec_G)))
    random.shuffle(permutation)
    k = Fr(random.randint(1, Fr.field_modulus))

    vec_R = [tracker.r_G for tracker in pre_shuffle_trackers]
    vec_S = [tracker.k_r_G for tracker in pre_shuffle_trackers]

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

    post_trackers = [WhiskTracker(r_G, k_r_G) for r_G, k_r_G in zip(vec_T, vec_U)]

    return post_trackers, shuffle_proof.to_json().encode()


SerializedWhiskTrackerProof = bytes


def IsValidWhiskOpeningProof(
    tracker: WhiskTracker,
    k_commitment: BLSG1Point,
    tracker_proof: SerializedWhiskTrackerProof,
) -> bool:
    """
    Verify knowledge of `k` such that `tracker.k_r_G == k * tracker.r_G` and `k_commitment == k * BLS_G1_GENERATOR`.
    """
    tracker_proof_instance = TrackerOpeningProof.from_json(tracker_proof.decode())

    transcript_verifier = CurdleproofsTranscript(b"whisk_opening_proof")
    return tracker_proof_instance.verify(
        transcript_verifier,
        affine_to_projective(tracker.k_r_G),
        affine_to_projective(tracker.r_G),
        affine_to_projective(k_commitment),
    )


def GenerateWhiskTrackerProof(
    tracker: WhiskTracker,
    k_G: BLSG1Point,
    k: Fr,
) -> SerializedWhiskTrackerProof:
    transcript_prover = CurdleproofsTranscript(b"whisk_opening_proof")
    opening_proof = TrackerOpeningProof.new(
        k_r_G=affine_to_projective(tracker.k_r_G),
        r_G=affine_to_projective(tracker.r_G),
        k_G=affine_to_projective(k_G),
        G=G1,
        k=k,
        transcript=transcript_prover,
    )

    return opening_proof.to_json().encode()
