import random
from typing import Container, Sequence, Tuple, NewType
from curdleproofs.crs import CurdleproofsCrs
from curdleproofs.curdleproofs import (
    CurdleProofsProof,
    shuffle_permute_and_commit_input,
)
from curdleproofs.curdleproofs_transcript import CurdleproofsTranscript
from curdleproofs.opening import TrackerOpeningProof
from curdleproofs.util import (
    PointProjective,
    affine_to_projective,
    Fr,
    BufReader,
)
from py_ecc.optimized_bls12_381.optimized_curve import G1, normalize, multiply
from py_ecc.bls.g2_primitives import (
    G1_to_pubkey,
    pubkey_to_G1,
)
from eth_typing import BLSPubkey

Bytes32 = NewType('Bytes32', bytes)

class WhiskTracker:
    r_G: BLSPubkey  # r * G
    k_r_G: BLSPubkey  # k * r * G

    def __init__(self, r_G: BLSPubkey, k_r_G: BLSPubkey):
        self.r_G = r_G
        self.k_r_G = k_r_G


SerializedCurdleProofsProof = bytes


def IsValidWhiskShuffleProof(
    crs: CurdleproofsCrs,
    pre_shuffle_trackers: Sequence[WhiskTracker],
    post_shuffle_trackers: Sequence[WhiskTracker],
    m: BLSPubkey,
    shuffle_proof: SerializedCurdleProofsProof,
) -> bool:
    """
    Verify `post_shuffle_trackers` is a permutation of `pre_shuffle_trackers`.
    """
    try:
        vec_R = [pubkey_to_G1(tracker.r_G) for tracker in pre_shuffle_trackers]
        vec_S = [pubkey_to_G1(tracker.k_r_G) for tracker in pre_shuffle_trackers]

        vec_T = [pubkey_to_G1(tracker.r_G) for tracker in post_shuffle_trackers]
        vec_U = [pubkey_to_G1(tracker.k_r_G) for tracker in post_shuffle_trackers]

        n = crs.n_el() + crs.n_blinders()
        shuffle_proof_instance = CurdleProofsProof.from_bytes(BufReader(shuffle_proof), n)
        M = pubkey_to_G1(m)

        shuffle_proof_instance.verify(crs, vec_R, vec_S, vec_T, vec_U, M)
        return True
    except:  # noqa: E722
        return False


def GenerateWhiskShuffleProof(
    crs: CurdleproofsCrs, pre_shuffle_trackers: Sequence[WhiskTracker]
) -> Tuple[Sequence[WhiskTracker], BLSPubkey, SerializedCurdleProofsProof]:
    permutation = list(range(crs.n_el()))
    random.shuffle(permutation)
    k = Fr(random.randint(1, Fr.field_modulus))

    vec_R = [pubkey_to_G1(tracker.r_G) for tracker in pre_shuffle_trackers]
    vec_S = [pubkey_to_G1(tracker.k_r_G) for tracker in pre_shuffle_trackers]

    vec_T, vec_U, M, vec_m_blinders = shuffle_permute_and_commit_input(
        crs, vec_R, vec_S, permutation, k
    )

    shuffle_proof = CurdleProofsProof.new(
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

    post_trackers = [WhiskTracker(G1_to_pubkey(r_G), G1_to_pubkey(k_r_G)) for r_G, k_r_G in zip(vec_T, vec_U)]

    return post_trackers, G1_to_pubkey(M), shuffle_proof.to_bytes()


SerializedWhiskTrackerProof = bytes


def IsValidWhiskOpeningProof(
    tracker: WhiskTracker,
    k_commitment: BLSPubkey,
    tracker_proof: SerializedWhiskTrackerProof,
) -> bool:
    """
    Verify knowledge of `k` such that `tracker.k_r_G == k * tracker.r_G` and `k_commitment == k * BLS_G1_GENERATOR`.
    """
    try:
        tracker_proof_instance = TrackerOpeningProof.from_bytes(BufReader(tracker_proof))

        transcript_verifier = CurdleproofsTranscript(b"whisk_opening_proof")
        tracker_proof_instance.verify(
            transcript_verifier,
            pubkey_to_G1(tracker.k_r_G),
            pubkey_to_G1(tracker.r_G),
            pubkey_to_G1(k_commitment),
        )
        return True
    except:  # noqa: E722
        return False


def GenerateWhiskTrackerProof(
    tracker: WhiskTracker,
    k: Fr,
) -> SerializedWhiskTrackerProof:
    transcript_prover = CurdleproofsTranscript(b"whisk_opening_proof")
    opening_proof = TrackerOpeningProof.new(
        k_r_G=pubkey_to_G1(tracker.k_r_G),
        r_G=pubkey_to_G1(tracker.r_G),
        k_G=multiply(G1, int(k)),
        k=k,
        transcript=transcript_prover,
    )

    return opening_proof.to_bytes()
