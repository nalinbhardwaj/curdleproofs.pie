from typing import Container, Sequence, Tuple
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.curdleproofs import N_BLINDERS, CurdleProofsProof
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from curdleproofs.opening import TrackerOpeningProof
from curdleproofs.util import PointAffine as BLSG1Point, affine_to_projective
from py_ecc.optimized_bls12_381.optimized_curve import G1


class WhiskTracker(Container):
    r_G: BLSG1Point  # r * G
    k_r_G: BLSG1Point  # k * r * G


SerializedCurdleProofsProof = str


def IsValidWhiskShuffleProof(
    pre_shuffle_trackers: Sequence[WhiskTracker],
    post_shuffle_trackers: Sequence[WhiskTracker],
    M: BLSG1Point,
    shuffle_proof: SerializedCurdleProofsProof,
) -> Tuple[bool, str]:
    """
    Verify `post_shuffle_trackers` is a permutation of `pre_shuffle_trackers`.
    """
    crs = CurdleproofsCrs.new(len(pre_shuffle_trackers), N_BLINDERS)
    vec_R = [pre_shuffle_tracker.r_G for pre_shuffle_tracker in pre_shuffle_trackers]
    vec_S = [pre_shuffle_tracker.k_r_G for pre_shuffle_tracker in pre_shuffle_trackers]

    vec_T = [post_shuffle_tracker.r_G for post_shuffle_tracker in post_shuffle_trackers]
    vec_U = [
        post_shuffle_tracker.k_r_G for post_shuffle_tracker in post_shuffle_trackers
    ]

    shuffle_proof_instance = CurdleProofsProof.from_json(shuffle_proof)

    M_projective = affine_to_projective(M)

    return shuffle_proof_instance.verify(crs, vec_R, vec_S, vec_T, vec_U, M_projective)


SerializedWhiskTrackerProof = str


def IsValidWhiskOpeningProof(
    tracker: WhiskTracker,
    k_commitment: BLSG1Point,
    tracker_proof: SerializedWhiskTrackerProof,
) -> bool:
    """
    Verify knowledge of `k` such that `tracker.k_r_G == k * tracker.r_G` and `k_commitment == k * BLS_G1_GENERATOR`.
    """
    tracker_proof_instance = TrackerOpeningProof.from_json(tracker_proof)

    if not tracker_proof_instance.G == G1:
        return False
    if not tracker_proof_instance.k_G == k_commitment:
        return False
    if not tracker_proof_instance.k_r_G == tracker.k_r_G:
        return False
    if not tracker_proof_instance.r_G == tracker.r_G:
        return False

    transcript_verifier = CurdleproofsTranscript(b"whisk_opening_proof")
    return tracker_proof_instance.verify(transcript_verifier)
